import ihm
import os

# If we're running from an SGE job, override the from_pubmed_id() function
# to return a cached value, since we don't have network access (needed to
# query PubMed directly)

def mock_from_pubmed(cls, pubmed_id):
    return ihm.Citation(
            pmid=27250206,
            title='Structure of Complement C3(H2O) Revealed By Quantitative '
                  'Cross-Linking/Mass Spectrometry And Modeling.',
            journal='Mol Cell Proteomics', volume=15, page_range=(2730,2743),
            year=2016, doi='10.1074/mcp.M115.056473', authors=[
                'Chen ZA', 'Pellarin R', 'Fischer L', 'Sali A', 'Nilges M',
                'Barlow PN', 'Rappsilber J'])

if 'JOB_ID' in os.environ:
    ihm.Citation.from_pubmed_id = classmethod(mock_from_pubmed)
