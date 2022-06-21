# SignalProcessingCup
Discrete candidate analysis for heart rate monitoring using wrist-type photoplethysmographic signals during intensive exercise

## Team
Kai Lin1, James Smith1, Ross Leon1, Robert Makepeace1, Nicholas Cummins1;2, Vidhyasaharan Sethu1
1: School of Elec. Eng. and Telecom., The University of New South Wales, Sydney Australia
2: ATP Research Lab, National ICT (NICTA), Australia

## Abstract
Wrist-type photoplethysmographic (PPG) signals are an increasingly popular way to monitor heart rate during intensive exercise. However these signals are highly influenced by motion artefacts. This paper proposes a novel system for heart rate extraction termed Discrete Candidate Analysis (DCA). The DCA extracts a discrete set of possible heart rates from PPG signal, from which the sequence of most likely candidates is chosen based on accelerometric data, harmonic and temporal analysis. On a dataset of 23 PPG recordings, the proposed system obtained an average error of 3.02 beats per minute.

![Brain CrossSection](/Assets/Brain.bmp)
Fig. 1. DCA Overall system: sets of heart rate and artefact candidates are chosen from PPG and accelerometer signals respectively. They are analysed and one heart rate candidate is chosen for output
## Report
[JournalArticle](/Assets/JournalArticle.pdf)

## Data
[Data](/Data)

## Software
[Software](/Software)

## Data
[Data](/Data)

## Images
[Results](/Results)

![Results](Results)
Fig. 2. Comparison of Average Absolute Error (AAE) between TROIKA [5], system of [4] and proposed DCA system over 23 datasets
