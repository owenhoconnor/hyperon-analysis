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


def comment_line(line: str) -> str:
    stripped = line.lstrip()
    if stripped.startswith("//"):
        return line
    indent = line[: len(line) - len(stripped)]
    return f"{indent}// {stripped}"


def remove_marked_block(text: str, start_marker: str, end_marker: str) -> str:
    pattern = re.compile(re.escape(start_marker) + r".*?" + re.escape(end_marker), re.DOTALL)
    return pattern.sub("", text)


def find_first(text, patterns):
    for pat in patterns:
        m = re.search(pat, text, flags=re.MULTILINE)
        if m:
            return m.group(1)
    return None


def comment_matching_lines(text: str, pattern: str) -> str:
    out = []
    regex = re.compile(pattern)
    for line in text.splitlines():
        if regex.search(line):
            out.append(comment_line(line))
        else:
            out.append(line)
    return "\n".join(out) + ("\n" if text.endswith("\n") else "")


def remove_matching_lines(text: str, pattern: str) -> str:
    return re.sub(pattern, "", text, flags=re.MULTILINE)


def replace_sig_legend_block(text: str, signif_base: str) -> str:
    lines = text.splitlines()
    sig_idx = None
    for i, line in enumerate(lines):
        if 'AddEntry(' in line and 'S/#sqrt{S+B}' in line and not line.lstrip().startswith("//"):
            sig_idx = i
            break

    new_block = [
        SIGLEG_START,
        '   leg = new TLegend(0.508, 0.8, 0.9, 0.92, nullptr, "brNDC");',
        '   leg->SetBorderSize(1);',
        '   leg->SetTextFont(62);',
        '   leg->SetTextSize(0.04);',
        '   leg->SetLineColor(TColor::GetColor("#7d8b9d"));',
        '   leg->SetLineStyle(1);',
        '   leg->SetLineWidth(1);',
        '   leg->SetFillColor(0);',
        '   leg->SetFillStyle(1);',
        '',
        f'   legentry = leg->AddEntry("{signif_base}","S/#sqrt{{S+B}}","L");',
        '   legentry->SetLineColor(TColor::GetColor("#00aa00"));',
        '   legentry->SetLineWidth(3);',
        '   legentry->SetTextFont(62);',
        '',
        '   leg->Draw("same");',
        SIGLEG_END,
    ]
    new_block_text = "\n".join(new_block)

    if sig_idx is None:
        m = re.search(r'^\s*TGaxis \*gaxis = new TGaxis\(', text, flags=re.MULTILINE)
        if not m:
            return text + "\n" + new_block_text + "\n"
        idx = m.start()
        prefix = "\n" if idx > 0 and text[idx - 1] != "\n" else ""
        return text[:idx] + prefix + new_block_text + "\n" + text[idx:]

    start = sig_idx
    while start >= 0 and 'leg = new TLegend(' not in lines[start]:
        start -= 1
    if start < 0:
        return text

    end = sig_idx
    while end < len(lines) and 'leg->Draw("same");' not in lines[end]:
        end += 1
    if end >= len(lines):
        return text

    replaced = lines[:start] + new_block + lines[end + 1:]
    return "\n".join(replaced) + ("\n" if text.endswith("\n") else "")


def ensure_axis_and_draw_block(text: str, sig_obj: str, bgd_obj: str, signif_obj: str, classifier: str) -> str:
    """
    Remove all existing sig/bgd/significance draw lines and all existing sig-axis lines,
    then inject one clean block that sets the axes on sigEffi and draws all three curves.
    """
    # Remove existing sig axis-setting lines and any draw lines (commented or uncommented)
    patterns = [
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetXaxis\(\)->SetTitle\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetXaxis\(\)->SetLabelOffset\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetXaxis\(\)->SetLabelSize\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetXaxis\(\)->SetTitleSize\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetXaxis\(\)->SetTitleOffset\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetXaxis\(\)->SetLabelFont\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetXaxis\(\)->SetTitleFont\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetYaxis\(\)->SetTitle\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetYaxis\(\)->SetLabelOffset\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetYaxis\(\)->SetLabelSize\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetYaxis\(\)->SetTitleSize\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetYaxis\(\)->SetTitleOffset\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetYaxis\(\)->SetLabelFont\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetYaxis\(\)->SetTitleFont\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetZaxis\(\)->SetLabelFont\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetZaxis\(\)->SetTitleOffset\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->GetZaxis\(\)->SetTitleFont\(.*?\);\n?',
        rf'^\s*//?\s*{re.escape(sig_obj)}->Draw\(".*?"\);\n?',
        rf'^\s*//?\s*{re.escape(bgd_obj)}->Draw\(".*?"\);\n?',
        rf'^\s*//?\s*{re.escape(signif_obj)}->Draw\(".*?"\);\n?',
    ]
    for pat in patterns:
        text = remove_matching_lines(text, pat)

    block = f'''   {sig_obj}->GetXaxis()->SetTitle("Cut value applied on {classifier} output");
   {sig_obj}->GetXaxis()->SetLabelOffset(0.012);
   {sig_obj}->GetXaxis()->SetLabelSize(0.03);
   {sig_obj}->GetXaxis()->SetTitleSize(0.04);
   {sig_obj}->GetXaxis()->SetTitleOffset(1.25);

   {sig_obj}->GetYaxis()->SetTitle("Efficiency");
   {sig_obj}->GetYaxis()->SetLabelOffset(0.01);
   {sig_obj}->GetYaxis()->SetLabelSize(0.03);
   {sig_obj}->GetYaxis()->SetTitleSize(0.04);
   {sig_obj}->GetYaxis()->SetTitleOffset(0.9);
   {sig_obj}->GetXaxis()->SetLabelFont(42);
   {sig_obj}->GetXaxis()->SetTitleOffset(1);
   {sig_obj}->GetXaxis()->SetTitleFont(42);
   {sig_obj}->GetYaxis()->SetLabelFont(42);
   {sig_obj}->GetYaxis()->SetTitleFont(42);
   {sig_obj}->GetZaxis()->SetLabelFont(42);
   {sig_obj}->GetZaxis()->SetTitleOffset(1);
   {sig_obj}->GetZaxis()->SetTitleFont(42);
   {sig_obj}->Draw("histl");
   {bgd_obj}->Draw("samehistl");
   {signif_obj}->Draw("samehistl");'''

    m = re.search(r'^\s*TLegend \*leg = new TLegend\(', text, flags=re.MULTILINE)
    if m:
        start = m.start()
        return text[:start] + block + "\n\n" + text[start:]

    return text + "\n" + block + "\n"


def patch_macro(path: Path):
    text = path.read_text(encoding="utf-8", errors="ignore")
    original = text

    text = remove_marked_block(text, SIGLEG_START, SIGLEG_END)
    text = remove_marked_block(text, SAVE_START, SAVE_END)

    classifier = CLASSIFIER_LABELS[path.name]
    output_base = OUTPUT_NAMES[path.name]
    var_num = path.parent.name.replace("-var", "")

    sig_obj = find_first(text, [
        r'TH1[FD]\s+\*(sigEffi__\d+)\s*=',
        r'(sigEffi__\d+)->Draw\(".*?"\);',
    ])
    bgd_obj = find_first(text, [
        r'TH1[FD]\s+\*(bgdEffi__\d+)\s*=',
        r'(bgdEffi__\d+)->Draw\(".*?"\);',
    ])
    signif_obj = find_first(text, [
        r'TH1[FD]\s+\*(significance_[A-Za-z0-9]+__\d+)\s*=',
        r'(significance_[A-Za-z0-9]+__\d+)->Draw\(".*?"\);',
    ])
    effpur_obj = find_first(text, [
        r'TH1[FD]\s+\*(effpurS_[A-Za-z0-9]+__\d+)\s*=',
        r'(effpurS_[A-Za-z0-9]+__\d+)->Draw\(".*?"\);',
    ])
    pur_obj = find_first(text, [
        r'TH1[FD]\s+\*(purS_[A-Za-z0-9]+__\d+)\s*=',
        r'(purS_[A-Za-z0-9]+__\d+)->Draw\(".*?"\);',
    ])
    sameaxis_obj = find_first(text, [
        r'(effpurS_[A-Za-z0-9]+__\d+)->Draw\("sameaxis"\);'
    ])
    canvas_name = find_first(text, [
        r'TCanvas \*(\w+)\s*=\s*new TCanvas\('
    ])

    if not (sig_obj and bgd_obj and signif_obj and canvas_name):
        print(f"Skipping {path}: could not identify sig/bgd/significance/canvas objects")
        return

    signif_base = re.sub(r'__\d+$', '', signif_obj)

    if effpur_obj:
        text = comment_matching_lines(text, rf'^\s*{re.escape(effpur_obj)}->Draw\(".*?"\);')
    if pur_obj:
        text = comment_matching_lines(text, rf'^\s*{re.escape(pur_obj)}->Draw\(".*?"\);')
    if sameaxis_obj:
        text = comment_matching_lines(text, rf'^\s*{re.escape(sameaxis_obj)}->Draw\("sameaxis"\);')

    # Comment ALL tex-related lines
    text = comment_matching_lines(text, r'^\s*TLatex\s+\*tex\s*=')
    text = comment_matching_lines(text, r'^\s*tex\s*=')
    text = comment_matching_lines(text, r'^\s*tex->')

    # Inject one clean axis+draw block for all three curves
    text = ensure_axis_and_draw_block(text, sig_obj, bgd_obj, signif_obj, classifier)

    # Replace significance legend cleanly
    text = replace_sig_legend_block(text, signif_base)

    # Match your example TGaxis styling
    text = re.sub(r'(gaxis->SetLabelSize\().*?(\);)', r'\g<1>0.03\2', text)
    text = re.sub(r'(gaxis->SetTitleSize\().*?(\);)', r'\g<1>0.04\2', text)
    text = re.sub(r'(gaxis->SetTitleOffset\().*?(\);)', r'\g<1>1\2', text)

    save_block = (
        f"{SAVE_START}\n"
        f'   {canvas_name}->SaveAs("{output_base}_{var_num}.pdf");\n'
        f"{SAVE_END}\n"
    )

    if SAVE_START not in text:
        text = re.sub(r'^\}', save_block + '}', text, flags=re.MULTILINE, count=1)

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
