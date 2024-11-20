# Consensus Update

[Consensus Script: Mult. Alleles By Length](llcombo6.py)

[Consensus Script: Mult. Alleles By Length Test File Output](llcombo6_output_nov_18.vcf)

[Consensus Script: Mult. Alleles By AD](llcombo9.py)

[Consensus Script: Mult. Alleles By AD Test File Output](llcombo9_output.vcf)

## Test Consensus File Stats (Consensus Script: Mult. Alleles By Length)
**Total calls written**: 7918


| Caller | Total Counts in Original VCF File | Total counts in Consensus | % Retained |
| --- | --- | --- | --- |
| LRSV | 4416 | 2843 | 64.39% |
| PBSV | 8435 | 2130 | 25.25% |
| Sniffles | 4906 | 2945 | 60.02% |


## Consensus Counts by Caller - Test Data
![alt text](https://github.com/Meshinchi-Lab/BGMP_Student_Project_2024/blob/f6b72304f270a7deec6828db82d7667026e71475/VCF_Consensus/testdata_consensus_bycaller.png)


## Consensus Counts by Caller - Full Data (llcombo8 output)
![consensus_counts](https://github.com/user-attachments/assets/fe4d1c42-007b-4270-98f0-0a32e9ef6e85)


## Full dataset pre-consensus stats
<img width="779"
                                                                                                                                          alt="Screenshot 2024-11-13 at 1 54 52 PM" src="https://github.com/user-attachments/assets/8ef43368-0c27-4055-90d5-86791f8510e7">

<img width="780" alt="Screenshot 2024-11-13 at 1 57 38 PM" src="https://github.com/user-attachments/assets/005dd77d-a721-4d37-bc8e-7c95f5976320">

## Assumptions/Rules of Our Algorithm
- Exclusively running on chromosomes 1-22
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

