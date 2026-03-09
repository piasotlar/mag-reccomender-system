from Bio import Entrez
from urllib.error import HTTPError, URLError
import pandas as pd
from collections import defaultdict
import ast
import time

Entrez.email = "ps02292@student.uni-lj.si"  
Entrez.api_key = "a5b4487a22d64bedc996512e7b495df47608"
GENOME_COVERAGE_EXCLUDE = "Partial genome (non-compliant)"
NCBI_EXCEL_PATH = "VMR_MSL40.v2.20260223.xlsx"
SEGMENTED = True


def nuccore_length(accession: str, retries: int = 4, base_delay: float = 0.3) -> int | None:
    accession = str(accession).strip()
    if not accession:
        return None

    for attempt in range(retries):
        try:
            h = Entrez.esearch(
                db="nuccore",
                term=f"{accession}[ACCN]",
                retmax=1
            )
            try:
                r = Entrez.read(h)
            finally:
                h.close()

            if not r.get("IdList"):
                return None

            uid = r["IdList"][0]

            h = Entrez.esummary(db="nuccore", id=uid)
            try:
                s = Entrez.read(h)
            finally:
                h.close()

            return int(s[0]["Length"])

        except (HTTPError, URLError, RuntimeError) as e:
            if attempt == retries - 1:
                raise
            time.sleep(base_delay * (2 ** attempt))


def collect_accessions(g):
    return g["Virus GENBANK accession"].dropna().tolist()

def ensure_list(x):
    if isinstance(x, list):
        return x
    if pd.isna(x):
        return []
    if isinstance(x, str):
        s = x.strip()
        # if it's a string that looks like a Python list: "['A','B']"
        if s.startswith("[") and s.endswith("]"):
            try:
                v = ast.literal_eval(s)
                return v if isinstance(v, list) else [v]
            except Exception:
                return [s]
        return [s]
    return [x]

def accessions_with_segment_dict(accessions_list):
    items = ensure_list(accessions_list)

    segmented = defaultdict(list)  # segment -> [accessions]

    for el in items:
        s = str(el).strip()

        if ";" not in s:
            if ":" in s:
                segment, acc = s.split(":", 1)
                seg = segment.strip()
                segmented[seg].append(acc)
            else:
                acc = s.split()[0] if "(" in s else s
                segmented[None].append(acc)
        else:
            for i, part in enumerate(s.split(";")):
                part = part.strip()
                if not part:
                    continue

                if ":" in part:
                    segment, acc = part.split(":", 1)
                    seg = segment.strip()
                    acc = acc.strip().strip(" ;,")
                    segmented[seg].append(acc)
                else:
                    acc = part.strip().strip(" ;,")
                    segmented[str(i+1)].append(acc)

    return dict(segmented)

def get_accessions(accessions_list):
    items = ensure_list(accessions_list)

    accessions = list()

    for el in items:
        s = str(el).strip()

        if ";" not in s:
            if ":" in s:
                _, acc = s.split(":", 1)
                acc = acc.strip()
            else:
                acc = s.split()[0] if "(" in s else s
            accessions.append(acc)
        else:
            for i, part in enumerate(s.split(";")):
                part = part.strip()
                if not part:
                    continue

                if ":" in part:
                    segment, acc = part.split(":", 1)
                    seg = segment.strip()
                    acc = acc.strip().strip(" ;,")
                    accessions.append(acc)
                else:
                    acc = part.strip().strip(" ;,")
                    accessions.append(acc)

    return accessions

def round_min_length(value):
    if value is not None and value != "":
        value = int(value)
        if value < 1000:
            value = (value // 100) * 100
        else:
            value = (value // 1000) * 1000

    else:
        value = ""
    return value

def round_max_length(value):
    if value is not None and value != "":
        value = int(value)
        if value < 1000:
            value = (value // 100) * 100
        else:
            value = (value // 1000) * 1000
    else:
        value = ""
    return value

def get_family_min_max(df_genuses):
    values = []
    values.extend(v for v in df_genuses["min_length_rounded"] if v != "")
    values.extend(v for v in df_genuses["max_length_rounded"] if v != "")
    if not values:
        return "", ""
    return min(values), max(values)



def main():
    df = pd.read_excel(NCBI_EXCEL_PATH)
    df = df[df["Genome coverage"] != GENOME_COVERAGE_EXCLUDE].copy()

    fg_accessions = (
        df.groupby(["Family", "Genus"], dropna=False)
        .apply(collect_accessions)
        .reset_index(name="accessions")
    )

    if SEGMENTED:
        fg_accessions["parsed_accessions"] = fg_accessions["accessions"].apply(accessions_with_segment_dict)
    else:
        fg_accessions["parsed_accessions"] = fg_accessions["accessions"].apply(get_accessions)

    failed = []
    rows = []

    for i, row in enumerate(fg_accessions.itertuples(index=False), start=1):
        family = row.Family
        genus = row.Genus
        accs = row.parsed_accessions

        if SEGMENTED:
            for segment, acc_list in accs.items():
                lengths = []

                for acc in acc_list:
                    try:
                        L = nuccore_length(acc) if acc else None
                    except (RuntimeError, HTTPError, URLError) as e:
                        failed.append((acc, str(e)))
                        L = None
                    except Exception as e:
                        failed.append((acc, f"{type(e).__name__}: {e}"))
                        L = None

                    if L is not None:
                        lengths.append(L)

                minimum = min(lengths) if lengths else ""
                maximum = max(lengths) if lengths else ""

                rows.append({
                    "Family": family,
                    "Genus": genus,
                    "Segment": segment,
                    "min_length": minimum,
                    "max_length": maximum,
                })
        else:
            lengths = []

            for acc in accs:
                try:
                    L = nuccore_length(acc) if acc else None
                except (RuntimeError, HTTPError, URLError) as e:
                    failed.append((acc, str(e)))
                    L = None
                except Exception as e:
                    failed.append((acc, f"{type(e).__name__}: {e}"))
                    L = None

                if L is not None:
                    lengths.append(L)

            minimum = min(lengths) if lengths else ""
            maximum = max(lengths) if lengths else ""

            rows.append({
                "Family": family,
                "Genus": genus,
                "min_length": minimum,
                "max_length": maximum,
            })

        if i % 10 == 0:
            print(f"processed {i}/{len(fg_accessions)} groups...")

    
    if failed:
        with open("failed_accessions.txt", "w", encoding="utf-8") as f:
            for acc, err in failed:
                f.write(f"{acc} -> {err}\n")

    fg_lengths = pd.DataFrame(rows)

    fg_lengths["min_length_rounded"] = fg_lengths["min_length"].apply(round_min_length)
    fg_lengths["max_length_rounded"] = fg_lengths["max_length"].apply(round_max_length)

    fg_lengths = fg_lengths.drop(columns=["min_length", "max_length"])

    if not SEGMENTED:
        family = (
            fg_lengths.groupby("Family")
            .apply(get_family_min_max)
            .reset_index(name="min_max")
        )

        family[["min_length_rounded", "max_length_rounded"]] = pd.DataFrame(
            family["min_max"].tolist(),
            index=family.index
        )
        family = family.drop(columns=["min_max"])

        host_map_family = (
            df.groupby("Family")["Host source"]
            .apply(lambda x: ", ".join(sorted(set(x.dropna().astype(str)))))
        )
        family["Host source"] = family["Family"].map(host_map_family)


    host_map_genus = (
        df.groupby("Genus")["Host source"]
        .apply(lambda x: ", ".join(sorted(set(x.dropna().astype(str)))))
    )
    
    fg_lengths["Host source"] = fg_lengths["Genus"].map(host_map_genus)
    

    if SEGMENTED:
        fg_lengths.to_csv(
            "genus_segmented_lengths_rounded_with_host.txt",
            sep="\t",
            index=False,
            na_rep="nan"
        )
    else:
        fg_lengths.to_csv(
            "genus_lengths_rounded_with_host.txt",
            sep="\t",
            index=False,
            na_rep="nan"
        )

        family.to_csv("family_lengths_rounded_with_host.txt", 
            sep="\t", 
            index=True, 
            na_rep="nan")

if __name__ == "__main__":
    main()