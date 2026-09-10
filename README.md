# The Role of CISD2 in RyR3-Mediated Calcium Release in Wolfram Syndrome Type 2

Bachelor's thesis project — Corso di laurea triennale in Biotecnologie Mediche, Università degli Studi di Milano (A.Y. 2024/2025), carried out at the **Laboratory of Molecular and Cellular Signaling, KU Leuven**.

- **Author:** Ludovico Salmin 
- **Supervisor:** Prof. Emanuela Galliera (UniMi)
- **Tutor:** Jens Loncke (KU Leuven)

## Background

Wolfram Syndrome is a rare neurodegenerative disorder with a life expectancy of around 30 years. Type 1 (WS1) is caused by mutations in *WFS1*, while Type 2 (WS2) is caused by mutations in *CISD2*. CISD2 is a transmembrane protein enriched at mitochondria-associated ER membranes (MAMs), containing a 2Fe-2S cluster domain, and is known to regulate Ca²⁺-mediated processes. Ryanodine receptors (RyRs) are the largest known ion channels, expressed on the ER of excitable tissues, with three isoforms (RyR1, RyR2, RyR3) activated by agonists such as ATP and caffeine.

## Project goals

- Investigate the physical interaction between CISD2 and RyR3
- Assess the impact of CISD2 on RyR3-mediated Ca²⁺ release
- Evaluate the effect of CISD2 mutants (CISD2^FeSCD, CISD2^N72S) on RyR3-mediated Ca²⁺ release
- Evaluate mitochondrial Ca²⁺ levels upon RyR3 stimulation

## Methods

- **Co-Immunoprecipitation:** C34 anti-RyR3 antibody, Dynabeads® Protein G, cell lysate samples, SDS-based elution buffer, to test the physical interaction between RyR3 and CISD2.
- **CRISPR-Cas9 genome editing:** HEK293 RyR3 cells transfected with pSpCas9(BB)-2A-Puro (PX459) carrying a CISD2-targeting gRNA, followed by puromycin selection, to generate CISD2 knockout cells.
- **Fluorescence microscopy:** cytosolic and mitochondrial Ca²⁺ measurements using FURA2-AM (excitation at 340/380 nm), mCherry (550 nm), and mitochondrial-targeted CEPIA (mtCEPIA, 470 nm).

## Results

- CISD2 co-precipitates with RyR3, indicating a physical interaction between the two proteins.
- CISD2 knockout cells show increased caffeine-induced RyR3-mediated Ca²⁺ release compared to control, suggesting CISD2 acts as a negative regulator of RyR3.
- CISD2^FeSCD and CISD2^N72S mutants do not reproduce the effect of full CISD2 loss, indicating these residues/domains are not solely responsible for the regulatory effect.
- In the absence of CISD2, mitochondrial Ca²⁺ uptake upon RyR3 stimulation is strongly reduced.

## Conclusions

- CISD2 is physically linked to RyR3.
- CISD2 acts as a negative regulator of RyR3-mediated Ca²⁺ release.
- The CISD2^FeSCD and CISD2^N72S mutants do not phenocopy the absence of CISD2.
- Without CISD2, Ca²⁺ fails to efficiently enter the mitochondria following RyR3 activation.

## Future directions

- Reverse Co-Immunoprecipitation to confirm the interaction from the CISD2 side
- Investigate the impact of CISD2 mutants on other components or pathways
- Increase the number of replicates for the mitochondrial Ca²⁺ studies

## Repository structure

```
.
├── plots/         # Generated figures and plots
├── presentation/  # Slides for the thesis defense
├── scripts/       # R scripts for statistic analysis and plots
└── thesis/        # Thesis manuscript
```

## Acknowledgments

ImageJ macros and Python scripts for processing of raw data were retrieved from https://github.com/jensloncke 
