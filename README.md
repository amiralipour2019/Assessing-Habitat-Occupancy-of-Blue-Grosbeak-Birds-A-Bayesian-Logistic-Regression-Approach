# Assessing Habitat Occupancy of Blue Grosbeak Birds: A Bayesian Logistic Regression Approach

## Introduction
This study focuses on the Blue Grosbeak birds' habitat occupancy in relation to field size. Given the potential for both observed and unobserved bird presence, we employ a Bayesian logistic regression model to analyze their habitat preferences. The goal is to provide insights into how field size impacts the likelihood of Blue Grosbeak occupancy, contributing valuable information for conservation efforts and habitat management.

## Dataset Description
The dataset comprises observations of Blue Grosbeak presence across various fields, encapsulated by the following features:
- **Field ID (s):** Unique identifier for each field location (s=1,2,...,41).
- **Field Size (X_s):** Numerical value representing the size of each field in the survey.
- **Occupancy Status (Y_{st}):** Binary indicator of the bird's occupancy status in field `s` during visit `t` (t=1,2,3), where 1 signifies presence and 0 absence.

## Methodology
The study adopts a Bayesian hierarchical model framework, leveraging the following assumptions:
- Bird occupancy status \(Y_{st}\) follows a Bernoulli distribution, influenced by a latent variable \(Z_s\) representing true occupancy.
- The probability of occupancy \(Z_s=1\) is modeled as a logistic function of field size, with parameters \(\alpha\) and \(\beta\) assumed to follow normal distributions with hyperparameter \(c\).

## Results and Analysis
We present the estimated parameters \(\alpha\), \(\beta\), and \(\theta\), highlighting their implications on habitat occupancy. The analysis includes:
- Interpretation of model parameters and their ecological significance.
- Evaluation of the model's primary assumptions and the posterior distributions of \(Z_s\) and \(\theta\).
- Implementation and findings from the Metropolis-Hasting algorithm to sample from the posterior distribution.

## Sensitivity Analysis of Hyperparameter \(c\)
We assess the impact of the hyperparameter \(c\) on model outcomes, aiming to identify its optimal value through comparative analysis.

## Conclusion
The study offers a nuanced understanding of Blue Grosbeak habitat preferences, demonstrating the utility of Bayesian logistic regression in ecological research. Our findings emphasize the significance of field size in determining habitat occupancy, with implications for conservation strategies.

## How to Run
Instructions on how to execute the analysis scripts, including dependencies, required packages, and execution commands.

## Acknowledgments
Acknowledgment of any contributions, data sources, or funding support.

## License
Details on the project's license and usage permissions.
