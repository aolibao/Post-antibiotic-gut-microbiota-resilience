#!/home/qqt/anaconda3/bin/python

from pathlib import Path
import re

ROOT = Path("/home/qqt/20061009response")
path = ROOT / "reviewer_response" / "Reviewer1_all_responses_revised_bilingual.md"
text = path.read_text(encoding="utf-8")

# Idempotent removal of locations inserted by an earlier run.
text = re.sub(
    r"\n\*\*Location in revised manuscript / 修订位置：\*\*[^\n]*\n",
    "\n",
    text,
)

locations = {
    "# Comment 1": "Lines 115–124 and 312–318; revised Figs. 2–4 / 第115–124行及312–318行；修订图2–4",
    "# Comment 2": "Lines 115–124 and 317–319 / 第115–124行及317–319行",
    "# Comment 3": "Lines 125–133 and 319–322 / 第125–133行及319–322行",
    "# Comment 4": "Lines 20–37 and 71–86; cohort flow figure / 第20–37行及71–86行；队列流程图",
    "# Comment 5": "Lines 71–86 and 410–413 / 第71–86行及410–413行",
    "# Major Comment 6": "Lines 134–151 and 175–183; revised longitudinal-analysis figure / 第134–151行及175–183行；修订纵向分析图",
    "# Major Comment 7": "Lines 142–151; revised Supplementary Fig. 3 legend / 第142–151行；修订补充图3图例",
    "# Major Comment 8": "Lines 98–109 and 357–364; revised Fig. 5 and Supplementary Fig. 4 terminology / 第98–109行及357–364行；修订图5及补充图4术语",
    "# Major Comment 9": "Lines 152–167, 324–336 and 365–384; revised Fig. 6 (lines 643–652) / 第152–167、324–336及365–384行；修订图6（第643–652行）",
    "# Major Comment 10": "Lines 2–3, 20–37 and 337–394 / 第2–3、20–37及337–394行",
    "## Minor Comment 1": "Lines 87–109; revised taxonomy labels and legends / 第87–109行；修订分类名称及图例",
    "## Minor Comment 2": "Lines 87–109 and 175–183 / 第87–109行及175–183行",
    "## Minor Comment 3": "Lines 134–141; revised Fig. 3 legend / 第134–141行；修订图3图例",
    "## Minor Comment 4": "Lines 134–141; revised Fig. 3 and companion supplementary heatmaps / 第134–141行；修订图3及补充线性热图",
    "## Minor Comment 5": "Revised figure legends at lines 580–652; vector PDF and 600-dpi TIFF files / 第580–652行图例；矢量PDF及600-dpi TIFF文件",
    "## Minor Comment 6": "Revised Supplementary Fig. 1 legend and axis labels / 修订补充图1图例及坐标轴",
    "## Minor Comment 7": "Lines 159–183 and revised boxplot legends / 第159–183行及修订箱线图图例",
    "## Minor Comment 8": "Lines 159–167 and 324–336; revised Supplementary Fig. 5 / 第159–167行及324–336行；修订补充图5",
    "## Minor Comment 9": "Lines 71–81 and 115–124 / 第71–81行及115–124行",
    "## Minor Comment 10": "Lines 410–413 / 第410–413行",
    "## Minor Comment 11": "Lines 20–37 / 第20–37行",
    "## Minor Comment 12": "Lines 82–86 and 152–158 / 第82–86行及152–158行",
}

ordered = list(locations)
for idx, heading in enumerate(ordered):
    start = text.find(heading + "\n")
    if start < 0:
        raise RuntimeError(f"Heading not found: {heading}")
    later = [text.find(h + "\n", start + len(heading) + 1) for h in ordered[idx + 1 :]]
    later = [x for x in later if x >= 0]
    end = min(later) if later else len(text)
    section = text[start:end]
    marker = "**中文回复：**"
    pos = section.find(marker)
    if pos < 0:
        raise RuntimeError(f"Chinese response marker not found in {heading}")
    insertion = (
        "\n**Location in revised manuscript / 修订位置：** "
        + locations[heading]
        + "\n\n"
    )
    section = section[:pos] + insertion + section[pos:]
    text = text[:start] + section + text[end:]

path.write_text(text, encoding="utf-8")
print(path)
