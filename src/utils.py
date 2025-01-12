# This file contains utility functions to perfomr analysis of the macrophytes.
# It includes tools to fetch ENTREZ, do stats, ...

from urllib.error import HTTPError
from Bio import Entrez
import altair as alt
import pandas as pd
import polars as pl


def get_taxon_id(taxon_name, ncbi_db="taxonomy", ncbi_idtype="acc"):
    """Function to obtain ncbi taxonomic ID of a given taxon name.

    Args:
        taxon_name (str): Query name of a taxon
        ncbi_db (str, optional): NCBI database that should be fetched. Defaults to "taxonomy",
        the taxonomic database of NCBI.
        ncbi_idtype (str, optional): Type of id that should be retrieved. Defaults to "acc", 
        the accession name.

    Returns:
        list: List of taxonomic ids of the taxon on NCBI.
    """
    try:
        with Entrez.esearch(db=ncbi_db, term=f'{taxon_name}[orgn]', retmax=10000) as handle:
            record = Entrez.read(handle)
    except HTTPError:
        return str('Database error, try later...')
    else:
        all_taxaIDs = record['IdList']

    return all_taxaIDs


def get_taxonlist_species(taxonlist, ncbi_db="taxonomy"):
    """Filter a list of taxa for their correspondence to be species or genus

    Args:
        taxonlist (list): List of taxon IDs (NCBI) that should be filtered.
        ncbi_db (str, optional): Database that should be fetched. Defaults to "taxonomy".

    Returns:
        List of scientific names of species and genera.
    """
    species_taxa = []
    for taxon_id in taxonlist:
        try:
            with Entrez.esummary(db=ncbi_db, id=taxon_id, retmax=10000) as handle:
                record = Entrez.read(handle)[0]
                
        except IndexError:
            print('New Error. Nemas Problemas!')

       # except HTTPError:
       #     print('New Error, but bro, stay seated, we skip it and keep searching...')

        else:
            if (record['Rank']=='species' and record['Genus'] != '') or record['Rank']=='genus':
                species_taxa.append(record["ScientificName"])

    return species_taxa
                
def plot_bar(df, category, sort_order, title, x_title, legend_title, filter_list=None):
    """Plotting categorical count data in a normalized barplot. 

    Args:
        df (polars dataframe): Dataframe that contains the data
        category (str): category that should be used for binning the data
        sort_order (list): order of entities of the category in which plot bars are plotted
        title (str): title of the figure
        x_title (str): title of the x axis
        legend_title (str): title of the legend
        filter_list (list, optional): List that filters the dataset for given entities of a category. Defaults to None.

    Returns:
        Barplot as a figure
    """
    df = df.group_by(category).len()
    df = df.to_pandas()
    df[category] = pd.Categorical(df[category], sort_order, ordered=True)
    df = pl.from_pandas(df)
    df = df.filter(pl.col(category).is_not_null())
    if filter_list is not None:
        df = df.filter(pl.col(category).is_in(filter_list))
    chart = df.plot.bar(
        x=alt.X("len", stack="normalize"),
        color=alt.Color(
            category, sort=alt.Sort(sort_order),
            legend=alt.Legend(title=legend_title)).scale(
                scheme="tableau20"
            ),
        order=alt.Order(f"color_{category}_sort_index:O"),

        ).properties(
        title=title
    )

    chart.encoding.x.title = f"{x_title}, N={df['len'].sum()}"
    return chart