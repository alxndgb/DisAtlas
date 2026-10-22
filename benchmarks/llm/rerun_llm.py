"""Opt-in re-run of the language-model comparators (used only by `score_llm.py --rerun-llm`).

Not needed to reproduce the paper: the default path scores the stored outputs in data/. This script re-sends
the shipped prompts with the shipped upload files and writes new DDA columns in the same layout as the stored
tables. A key is needed only here: set ANTHROPIC_API_KEY (Claude) and/or GEMINI_API_KEY (Gemini, and set its
model ID in models.json). Model outputs are not deterministic, so re-runs will not reproduce the paper's
numbers exactly. The paper's runs were made in each system's own interface (file upload + prompt, fresh session);
here the file is sent in chunks as tab-separated text with the same prompt, which is an approximation.

    python benchmarks/llm/score_llm.py --rerun-llm      # re-run, then score the new runs
    python benchmarks/llm/rerun_llm.py                  # re-run only
"""
import io
import json
import os
import re

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "rerun")
CFG = json.load(open(os.path.join(HERE, "models.json")))
FORMAT_NOTE = ("\n\nReturn the rows as tab-separated text with a header, with the same columns as given plus a final "
               "DDA column. Return only the table.")


def call_model(provider, model, prompt, max_tokens):
    if provider == "anthropic":
        import anthropic  # pip install anthropic
        r = anthropic.Anthropic().messages.create(model=model, max_tokens=max_tokens,
                                                  messages=[{"role": "user", "content": prompt}])
        return "".join(b.text for b in r.content if getattr(b, "type", "") == "text")
    if provider == "google":
        from google import genai  # pip install google-genai
        if not model:
            raise SystemExit("set systems.gemini_pro_3_1.model in models.json")
        return genai.Client().models.generate_content(model=model, contents=prompt).text
    raise ValueError(provider)


def parse_dda(text, ids):
    """Parse a returned table; DDA per blind_id, clipped to [0, 1]; missing rows stay NaN."""
    body = re.sub(r"^```[a-z]*\n|```\s*$", "", text.strip(), flags=re.M)
    try:
        t = pd.read_csv(io.StringIO(body), sep="\t")
        t = t.set_index("blind_id")["DDA"]
    except Exception:
        t = pd.Series({m.group(1): float(m.group(2)) for m in
                       re.finditer(r"(IXB\d+)\b.*?([01](?:\.\d+)?)\s*$", body, flags=re.M)})
    v = pd.to_numeric(t, errors="coerce").clip(0, 1)
    return v.reindex(ids).values


def run_system(sys_key, upload, prompt_file):
    s = CFG["systems"][sys_key]; req = CFG["request"]
    up = pd.read_csv(os.path.join(HERE, upload), sep="\t")
    prompt = open(os.path.join(HERE, prompt_file)).read()
    out = []
    for i in range(0, len(up), req["rows_per_request"]):
        chunk = up.iloc[i:i + req["rows_per_request"]]
        msg = prompt + FORMAT_NOTE + "\n\nAttached file:\n" + chunk.to_csv(sep="\t", index=False)
        out.append(parse_dda(call_model(s["provider"], s["model"], msg, req["max_tokens"]), chunk.blind_id))
    return up.blind_id, np.concatenate(out)


def rerun_all():
    os.makedirs(OUT, exist_ok=True)
    p2 = CFG["protocol2"]; ids = None; runs = {}
    for k, s in CFG["systems"].items():
        for r in range(1, s["runs"] + 1):
            ids, runs[f"{k}_run{r}"] = run_system(k, p2["upload"], p2["prompt"])
    R = pd.DataFrame({"blind_id": ids, **runs})
    paths = {"runs": os.path.join(OUT, "llm_run_scores.tsv.gz")}
    R.to_csv(paths["runs"], sep="\t", index=False)
    return paths


if __name__ == "__main__":
    print(rerun_all())
