# <a href="https://drugsearch.org"><img src="frontend/src/assets/drugsearch_logo.png" alt="drug_search_logo" width="200px" height="auto" /></a>

**DrugSearch** is a comprehensive tool designed to visualize the adverse effects of various drugs. Not only can users explore the potential side effects of a specific drug, but they can also work the other way around: inputting an adverse effect to discover which drugs might cause it.

## Features

- **Drug to Adverse Effect**: Enter a drug name to view its potential adverse effects.
- **Adverse Effect to Drug**: Discover which drugs can cause a specific adverse effect.
- **Visual Insights**: Utilize our visualizing tools powered by ApexCharts to get a graphical representation of the data.
- **Drug Information**: Fetch detailed information about a drug, including its molecule picture, IUPAC name, drug class, and more.
- **AI-Powered Summaries**: Generate a concise AI summary for quick statistics on any drug.
- **Reliable Data Sources**: Our data is fetched from the FDA API and a local database derived from DrugBank.ca, ensuring accuracy and comprehensiveness.

## Technologies Used

- **Frontend**: React
- **Visualization**: ApexCharts
- **Backend**: Flask
- **Data Sources**: FDA API and DrugBank.ca

## Usage

Access DrugSearch [here](https://drugsearch.org).

1. **Search**: Use the search bar to enter a drug name or an adverse effect.
2. **Drug Details**: Take a look at the chemical information or use the accordion next to it to unveil the pharmacological data.
3. **Summary**: Click on a book icon next to the drug name to reveal the summary.
4. **Visualize**: Scroll down to see the graphs.

## Acknowledgments

- Special thanks to [DrugBank.ca](https://www.drugbank.ca/) for providing a comprehensive database on drugs.
- [OpenFDA API](https://open.fda.gov/) for their reliable and extensive data.
