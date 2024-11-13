# Consensus Update

[Script](https://github.com/Meshinchi-Lab/BGMP_Student_Project_2024/blob/d0a516e8edc80ce6c213546a8bd3b27c9d9e1884/VCF_Consensus/llcombo5.py)

[Output File ](https://github.com/Meshinchi-Lab/BGMP_Student_Project_2024/blob/d0a516e8edc80ce6c213546a8bd3b27c9d9e1884/VCF_Consensus/testconsensus_llcombo5.vcf)

## Test Consensus File Stats
**Total calls written**: 7918


| Caller | Total Counts in Original VCF File | Total counts in Consensus | % Retained |
| --- | --- | --- | --- |
| LRSV | 4416 | 2843 | 64.39% |
| PBSV | 8435 | 2130 | 25.25% |
| Sniffles | 4906 | 2945 | 60.02% |


![alt text](https://github.com/Meshinchi-Lab/BGMP_Student_Project_2024/blob/f6b72304f270a7deec6828db82d7667026e71475/VCF_Consensus/testdata_consensus_bycaller.png)


## Assumptions/Rules of Our Algorithm
- If two or more calls from different callers overlap, the one with the highest allelic depth (AD) is written to the output file
- If the two or more calls from different callers overlap and have the same allelic depth, the first encountered call is written to the output file (“caller agnostic”)
- **Handling breakends**:
  - Established nomenclature for BND “type” per [VCF 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf):
    - BND_type = 1: N[
    - BND_type = 2: ]N
    - BND_type = 3: N]
    - BND_type = 4: [N
  - Considered a match between callers if:
    - BNDs are of the same type and have the same mate chromosome
    - BNDs are of type 2 and 4 and mate position of the type 2 breakend is >= the mate position of the type 4 breakend

