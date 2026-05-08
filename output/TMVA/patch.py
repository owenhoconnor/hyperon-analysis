from pathlib import Path
import re

FOLDERS = ["5-var", "10-var", "20-var", "32-var", "44-var"]
MACROS = ["rulefit.C", "bdt.C", "svm.C", "mlp.C", "knn.C", "ld.C", "pca.C", "likelihood.C"]

SIGLEG_START = "// BEGIN_AUTO_SIG_LEGEND"
SIGLEG_END = "// END_AUTO_SIG_LEGEND"
SAVE_START = "// BEGIN_AUTO_SAVEAS"
SAVE_END = "// END_AUTO_SAVEAS"

CLASSIFIER_LABELS = {
    "rulefit.C": "RuleFit",
    "bdt.C": "BDT",
    "svm.C": "SVM",
    "mlp.C": "MLPBNN",
    "knn.C": "KNN",
    "ld.C": "LD",
    "pca.C": "LikelihoodPCA",
    "likelihood.C": "Likelihood",
}

OUTPUT_NAMES = {
    "rulefit.C": "effs_RuleFit",
    "bdt.C": "effs_BDT",
    "svm.C": "effs_SVM",
    "mlp.C": "effs_MLPBNN",
    "knn.C": "effs_KNN",
    "ld.C": "effs_LD",
    "pca.C": "effs_LikelihoodPCA",
    "likelihood.C": "effs_Likelihood",
}

def comment_matching_lines(text: str, pattern: str) -> str:
    out = []
    regex = re.compile(pattern)
    for line in text.splitlines():
        stripped = line.lstrip()
        if stripped.startswith("//"):
            out.append(line)
            continue
        if regex.search(line):
            indent = line[: len(line) - len(line.lstrip())]
            out.append(f"{indent}// {line.lstrip()}")
        else:
            out.append(line)
    return "\n".join(out) + ("\n" if text.endswith("\n") else "")

def replace_or_insert_block(text, start_marker, end_marker, new_block, anchor_pattern=None, insert_after=False):
    block_pattern = re.compile(
        re.escape(start_marker) + r".*?" + re.escape(end_marker),
        re.DOTALL
    )
    if block_pattern.search(text):
        return block_pattern.sub(new_block, text, count=1)

    if anchor_pattern is None:
        return text + "\n" + new_block + "\n"

    anchor = re.search(anchor_pattern, text, flags=re.MULTILINE)
    if not anchor:
        return text

    idx = anchor.end() if insert_after else anchor.start()
    prefix = "\n" if idx > 0 and text[idx - 1] != "\n" else ""
    return text[:idx] + prefix + new_block + "\n" + text[idx:]

def find_first(text, patterns):
    for pat in patterns:
        m = re.search(pat, text, flags=re.MULTILINE)
        if m:
            return m.group(1)
    return None

def comment_extra_legends(text: str) -> str:
    """
    Keep the first live legend block (main legend), comment out later live legend blocks.
    """
    lines = text.splitlines()
    blocks = []
    i = 0
    while i < len(lines):
        if "TLegend *leg = new TLegend(" in lines[i] and not lines[i].lstrip().startswith("//"):
            start = i
            j = i
            while j < len(lines) and 'leg->Draw("same");' not in lines[j]:
                j += 1
            if j < len(lines):
                blocks.append((start, j))
                i = j + 1
                continue
        i += 1

    if len(blocks) <= 1:
        return text

    out = lines[:]
    for start, end in blocks[1:]:
        for k in range(start, end + 1):
            stripped = out[k].lstrip()
            if not stripped.startswith("//"):
                indent = out[k][: len(out[k]) - len(out[k].lstrip())]
                out[k] = f"{indent}// {out[k].lstrip()}"

    return "\n".join(out) + ("\n" if text.endswith("\n") else "")

def patch_macro(path: Path):
    text = path.read_text(encoding="utf-8", errors="ignore")
    original = text

    classifier = CLASSIFIER_LABELS[path.name]
    output_base = OUTPUT_NAMES[path.name]
    var_folder = path.parent.name  # e.g. 5-var
    var_num = var_folder.replace("-var", "")

    sig_obj = find_first(text, [
        r'(sigEffi__\d+)->Draw\(".*?"\);',
        r'TH1[FD]\s+\*(sigEffi__\d+)\s*=',
    ])
    bgd_obj = find_first(text, [
        r'(bgdEffi__\d+)->Draw\(".*?"\);',
        r'TH1[FD]\s+\*(bgdEffi__\d+)\s*=',
    ])
    signif_obj = find_first(text, [
        r'(significance_[A-Za-z0-9]+__\d+)->Draw\(".*?"\);',
        r'TH1[FD]\s+\*(significance_[A-Za-z0-9]+__\d+)\s*=',
    ])
    effpur_obj = find_first(text, [
        r'(effpurS_[A-Za-z0-9]+__\d+)->Draw\(".*?"\);',
        r'TH1[FD]\s+\*(effpurS_[A-Za-z0-9]+__\d+)\s*=',
    ])
    pur_obj = find_first(text, [
        r'(purS_[A-Za-z0-9]+__\d+)->Draw\(".*?"\);',
        r'TH1[FD]\s+\*(purS_[A-Za-z0-9]+__\d+)\s*=',
    ])
    sameaxis_obj = find_first(text, [
        r'(effpurS_[A-Za-z0-9]+__\d+)->Draw\("sameaxis"\);',
        r'TH1[FD]\s+\*(effpurS_[A-Za-z0-9]+__\d+)\s*=',
    ])
    canvas_name = find_first(text, [
        r'TCanvas \*(\w+)\s*=\s*new TCanvas\('
    ])

    if not (sig_obj and bgd_obj and signif_obj and canvas_name):
        print(f"Skipping {path}: could not identify sig/bgd/significance/canvas objects")
        return

    signif_base = re.sub(r'__\d+$', '', signif_obj)

    # 1) Comment out unwanted curves
    if effpur_obj:
        text = comment_matching_lines(text, rf'^{re.escape(effpur_obj)}->Draw\(".*?"\);')
    if pur_obj:
        text = comment_matching_lines(text, rf'^{re.escape(pur_obj)}->Draw\(".*?"\);')
    if sameaxis_obj:
        text = comment_matching_lines(text, rf'^{re.escape(sameaxis_obj)}->Draw\("sameaxis"\);')

    # 2) Comment out TLatex summary text fully
    text = comment_matching_lines(text, r'^\s*TLatex\s+\*tex\s*=\s*new TLatex\(')
    text = comment_matching_lines(text, r'^\s*tex\s*=\s*new TLatex\(')
    text = comment_matching_lines(text, r'^\s*tex->SetNDC\(\);')
    text = comment_matching_lines(text, r'^\s*tex->SetTextSize\(')
    text = comment_matching_lines(text, r'^\s*tex->Draw\(\);')

    # 3) Normalize live draw calls
    text = re.sub(
        rf'^\s*{re.escape(sig_obj)}->Draw\(".*?"\);',
        f'   {sig_obj}->Draw("histl");',
        text,
        flags=re.MULTILINE
    )
    text = re.sub(
        rf'^\s*{re.escape(bgd_obj)}->Draw\(".*?"\);',
        f'   {bgd_obj}->Draw("samehistl");',
        text,
        flags=re.MULTILINE
    )
    text = re.sub(
        rf'^\s*{re.escape(signif_obj)}->Draw\(".*?"\);',
        f'   {signif_obj}->Draw("samehistl");',
        text,
        flags=re.MULTILINE
    )

    # 4) Force x/y axis styling to match your example rulefit.C exactly
    replacements = [
        (rf'({re.escape(sig_obj)}->GetXaxis\(\)->SetTitle\().*?(\);)',
         rf'\1"Cut value applied on {classifier} output"\2'),
        (rf'({re.escape(sig_obj)}->GetXaxis\(\)->SetLabelOffset\().*?(\);)',
         rf'\g<1>0.012\2'),
        (rf'({re.escape(sig_obj)}->GetXaxis\(\)->SetLabelSize\().*?(\);)',
         rf'\g<1>0.03\2'),
        (rf'({re.escape(sig_obj)}->GetXaxis\(\)->SetTitleSize\().*?(\);)',
         rf'\g<1>0.04\2'),
        (rf'({re.escape(sig_obj)}->GetXaxis\(\)->SetTitleOffset\().*?(\);)',
         rf'\g<1>1.25\2'),
        (rf'({re.escape(sig_obj)}->GetYaxis\(\)->SetTitle\().*?(\);)',
         rf'\1"Efficiency"\2'),
        (rf'({re.escape(sig_obj)}->GetYaxis\(\)->SetLabelOffset\().*?(\);)',
         rf'\g<1>0.01\2'),
        (rf'({re.escape(sig_obj)}->GetYaxis\(\)->SetLabelSize\().*?(\);)',
         rf'\g<1>0.03\2'),
        (rf'({re.escape(sig_obj)}->GetYaxis\(\)->SetTitleSize\().*?(\);)',
         rf'\g<1>0.04\2'),
        (rf'({re.escape(sig_obj)}->GetYaxis\(\)->SetTitleOffset\().*?(\);)',
         rf'\g<1>0.9\2'),
        (rf'({re.escape(sig_obj)}->GetXaxis\(\)->SetLabelFont\().*?(\);)',
         rf'\g<1>42\2'),
        (rf'({re.escape(sig_obj)}->GetXaxis\(\)->SetTitleFont\().*?(\);)',
         rf'\g<1>42\2'),
        (rf'({re.escape(sig_obj)}->GetYaxis\(\)->SetLabelFont\().*?(\);)',
         rf'\g<1>42\2'),
        (rf'({re.escape(sig_obj)}->GetYaxis\(\)->SetTitleFont\().*?(\);)',
         rf'\g<1>42\2'),
    ]

    for pat, repl in replacements:
        text = re.sub(pat, repl, text)

    # 5) Keep only first live legend block, then add one significance legend
    text = comment_extra_legends(text)

    siglegend_block = f"""{SIGLEG_START}
   leg = new TLegend(0.508, 0.8, 0.9, 0.92, nullptr, "brNDC");
   leg->SetBorderSize(1);
   leg->SetTextFont(62);
   leg->SetTextSize(0.04);
   leg->SetLineColor(TColor::GetColor("#7d8b9d"));
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1);

   legentry = leg->AddEntry("{signif_base}","S/#sqrt{{S+B}}","L");
   legentry->SetLineColor(TColor::GetColor("#00aa00"));
   legentry->SetLineWidth(3);
   legentry->SetTextFont(62);

   leg->Draw("same");
{SIGLEG_END}"""

    text = replace_or_insert_block(
        text,
        SIGLEG_START,
        SIGLEG_END,
        siglegend_block,
        anchor_pattern=r'^\s*TGaxis \*gaxis = new TGaxis\(',
        insert_after=False,
    )

    # 6) Match TGaxis formatting from your example
    text = re.sub(r'(gaxis->SetLabelSize\().*?(\);)', r'\g<1>0.03\2', text)
    text = re.sub(r'(gaxis->SetTitleSize\().*?(\);)', r'\g<1>0.04\2', text)
    text = re.sub(r'(gaxis->SetTitleOffset\().*?(\);)', r'\g<1>1\2', text)

    # 7) SaveAs block
    save_block = f"""{SAVE_START}
   {canvas_name}->SaveAs("{output_base}_{var_num}.pdf");
{SAVE_END}"""

    text = replace_or_insert_block(
        text,
        SAVE_START,
        SAVE_END,
        save_block,
        anchor_pattern=r'^\}',
        insert_after=False,
    )

    if text != original:
        path.write_text(text, encoding="utf-8")
        print(f"Patched: {path}")
    else:
        print(f"No changes made: {path}")

def main():
    for folder in FOLDERS:
        folder_path = Path(folder)
        if not folder_path.exists():
            print(f"Folder not found: {folder}")
            continue

        for macro in MACROS:
            path = folder_path / macro
            if path.exists():
                patch_macro(path)
            else:
                print(f"Missing file: {path}")

if __name__ == "__main__":
    main()
