"""Optional: re-ask the model with the shipped prompts. Not used by any rebuild script.

Model output is not deterministic, so a re-run will not reproduce the frozen tables row for row
(re-judging 944 pairs gave 800/944 keep/drop agreement). The rebuild reads the frozen tables in
data/00_raw/llm_annotations/.

Opt-in only (the `--rerun-llm` path): needs the `anthropic` package and your own ANTHROPIC_API_KEY.
Reproducing the paper needs neither. Inputs are user-supplied CSVs, because disease names and drug
lists are not shipped:

    python call_llm.py judge  pairs.csv    out.jsonl   # columns as in prompts/ixidn_judge_user.txt
    python call_llm.py labels diseases.csv out.jsonl   # columns as in prompts/disease_label_user.txt

The label run in the paper also gave the model literature-search, ontology-lookup and URL-fetch
tools; this script offers only the submit_labels tool, so the model answers from its own knowledge.
"""
import json
import sys
from pathlib import Path

import anthropic
import pandas as pd

HERE = Path(__file__).resolve().parent
M = json.loads((HERE / "models.json").read_text())


def judge(df, out):
    cfg = M["ixidn_judge"]
    system = (HERE / cfg["system_prompt"]).read_text()
    line = (HERE / cfg["user_prompt"]).read_text().splitlines()[1]
    client = anthropic.Anthropic()
    rows = df.to_dict("records")
    with open(out, "w") as f:
        for i in range(0, len(rows), 20):
            msg = "Pairs:\n" + "\n".join(line.format(**r) for r in rows[i:i + 20])
            r = client.messages.create(model=cfg["model"], max_tokens=16000, system=system,
                                       messages=[{"role": "user", "content": msg}])
            f.write(json.dumps({"text": r.content[0].text}) + "\n")


def labels(df, out):
    cfg = M["disease_labels"]
    system = (HERE / cfg["system_prompt"]).read_text()
    tmpl = (HERE / cfg["user_prompt"]).read_text()
    tool = json.loads((HERE / cfg["tool"]).read_text())
    client = anthropic.Anthropic()
    with open(out, "w") as f:
        for r in df.to_dict("records"):
            res = client.messages.create(
                model=cfg["model"], max_tokens=4000, system=system, tools=[tool],
                tool_choice={"type": "tool", "name": tool["name"]},
                messages=[{"role": "user", "content": tmpl.format(**r)}])
            ans = next(b.input for b in res.content if b.type == "tool_use")
            f.write(json.dumps({"disease_id": r["disease_id"], **ans}) + "\n")


if __name__ == "__main__":
    task, inp, out = sys.argv[1:4]
    {"judge": judge, "labels": labels}[task](pd.read_csv(inp), out)
