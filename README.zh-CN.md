# Molecular Cloning Validation Skills 中文说明

这是一个用于分子克隆验证的可复用 AI agent skill 仓库。仓库里包含两个版本：一个给 Codex 使用，一个给 Claude Code 使用；两个版本共用同一套 Python pipeline 思路。

- `codex/molecular-cloning-validation/`: Codex 版本。
- `claude-code/molecular-cloning-validation/`: Claude Code 版本。

这个 pipeline 可以从合成序列表、受体载体图谱和可选的 Sanger 测序结果出发，生成 insert 图谱、克隆后质粒图谱、测序比对表格、Markdown 报告和 HTML 可视化报告。

## 具体流程

1. 读取 `.xlsx`、`.csv` 或 `.tsv` 格式的 construct 表格。
2. 读取受体载体文件，支持 SnapGene `.dna`、GenBank 或 FASTA。
3. 定位需要被替换的载体区域。可以用两种方式：
   - 用 `--replace-start` 和 `--replace-end` 提供准确的 1-based 坐标；
   - 用 `--left-anchor` 和 `--right-anchor` 提供替换区域左右边界序列。
4. 对表格中每条合成 insert，生成一个新的克隆后质粒序列。
5. 输出合成 insert 的 GenBank 文件和克隆后完整质粒 GenBank 文件。
6. 如果本机安装了 SnapGene CLI，则把 GenBank 转换成 SnapGene 原生 `.dna` 文件。
7. 如果有第二个固定 Gibson 片段，例如 tag、reporter、binding module 或 linker，则从原始载体 feature 中直接提取参考序列。
8. 将 Sanger read 比对到 `variable insert + fixed module` 的联合参考序列。
9. 输出每条 read 的差异、每个 construct 的 pass/fail 汇总和推荐 colony/sample。

## Anchor 和 Overlap 参数不是固定的

这些参数每个用户、每个载体、每个克隆设计都可能不一样：

- `--left-anchor`: 原始载体中替换区域左边界的序列。
- `--right-anchor`: 原始载体中替换区域右边界的序列。
- `--left-overlap-seq`: 合成 insert 中包含的 5' 同源臂或 overlap 序列。
- `--right-overlap-seq`: 合成 insert 中包含的 3' 同源臂、linker 或右侧边界序列。

所以你前面看到的具体序列只是某一个项目的例子，不是这个 skill 固定要求的序列。后续换成别人的载体或别人的 Gibson 设计时，agent 应该先读取用户提供的载体和合成序列，再确定这些参数。

如果已经知道准确替换坐标，建议直接用坐标，而不是 anchor：

```bash
--replace-start START_POS --replace-end END_POS
```

这里的坐标是 1-based、闭区间。

## 常用命令模板

```bash
python scripts/cloning_validation_pipeline.py \
  --table synthesis.xlsx \
  --name-column construct_id \
  --dna-column DNA_final \
  --name-map A:targetA,B:targetB \
  --vector recipient_vector.dna \
  --left-anchor LEFT_VECTOR_BOUNDARY_SEQUENCE \
  --right-anchor RIGHT_VECTOR_BOUNDARY_SEQUENCE \
  --left-overlap-seq INSERT_5P_OVERLAP_SEQUENCE \
  --right-overlap-seq INSERT_3P_OVERLAP_SEQUENCE \
  --second-feature-name "fixed module" \
  --sanger-success-dir sanger/success \
  --sanger-failed-dir sanger/failed \
  --out clone_analysis
```

如果没有安装 SnapGene CLI，可以加：

```bash
--skip-snapgene
```

这样仍然会生成 GenBank 文件和 Sanger 报告，只是不生成 `.dna` 文件。

## macOS 和 Windows 的 SnapGene CLI 是否一致

SnapGene CLI 是 SnapGene 官方桌面软件自带的命令行程序。macOS 和 Windows 上的核心参数基本一致，主要差别是可执行文件路径和 shell 换行/引号写法。

macOS 常见路径：

```bash
/Applications/SnapGene.app/Contents/MacOS/SnapGene
```

Windows 常见路径：

```powershell
C:\Program Files\SnapGene\SnapGene.exe
C:\Program Files (x86)\SnapGene\SnapGene.exe
```

脚本会自动尝试查找 `SnapGene`、`SnapGene.exe` 以及这些常见路径。如果自动查找失败，可以手动指定：

```bash
--snapgene-cli /Applications/SnapGene.app/Contents/MacOS/SnapGene
```

Windows PowerShell 示例：

```powershell
python scripts/cloning_validation_pipeline.py `
  --table synthesis.xlsx `
  --name-column construct_id `
  --dna-column DNA_final `
  --vector recipient_vector.dna `
  --replace-start START_POS `
  --replace-end END_POS `
  --second-feature-name "fixed module" `
  --snapgene-cli "C:\Program Files\SnapGene\SnapGene.exe" `
  --out clone_analysis
```

需要注意：某些 SnapGene 版本可能已经写出了有效 `.dna` 文件，但进程退出时仍然有 warning 或非零退出码。因此最好用 `snapgene_reader` 把生成的 `.dna` 文件读回，确认序列和预期一致。

## 输出目录

pipeline 会生成：

- `01_synthesis_genbank/`: 合成 insert 的 GenBank 文件。
- `02_synthesis_snapgene_dna/`: 合成 insert 的 SnapGene `.dna` 文件，如果转换成功。
- `03_cloned_plasmid_genbank/`: 克隆后完整质粒 GenBank 文件。
- `04_cloned_plasmid_snapgene_dna/`: 克隆后完整质粒 SnapGene `.dna` 文件，如果转换成功。
- `05_sanger_alignment/sanger_results.csv`: 每条 Sanger read 的详细差异。
- `05_sanger_alignment/sanger_summary.csv`: 每个 construct 的汇总和推荐样品。
- `05_sanger_alignment/sanger_report.md`: Markdown 报告。
- `05_sanger_alignment/sanger_report.html`: HTML 可视化报告。
- `05_sanger_alignment/second_fragment_reference.fasta`: 从载体 feature 中提取的固定片段参考序列。
- `generation_summary.csv`: insert 长度和克隆后质粒长度汇总。

## 依赖

```bash
pip install pandas openpyxl biopython snapgene-reader
```

如果需要生成 SnapGene 原生 `.dna` 文件，还需要本机安装 SnapGene 桌面软件。

## 安装为 Skill

Codex 使用时，把这个目录复制到 Codex skills 目录：

```text
codex/molecular-cloning-validation
```

Claude Code 使用时，把这个目录复制到 Claude Code skills 目录：

```text
claude-code/molecular-cloning-validation
```

## 隐私说明

这个仓库只包含可复用 workflow、skill 文档和脚本，不包含任何项目特异的序列、Sanger 原始数据或分析结果。
