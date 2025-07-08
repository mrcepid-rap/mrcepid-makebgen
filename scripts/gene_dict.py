"""
This script is used to create a gene dictionary from the Ensembl Biomart database, which contains information about genes, transcripts, and their genomic locations.
This is necessary for running the makebgen pipeline, as it allows for efficient chunking of genomic data based on gene locations.

Note that you need to generate this dictionary only once, and it should be stored on DNA Nexus (so that you can supply the DXid of the gene dictionary to the makebgen applet).
"""

import json

import pandas as pd
from pybiomart import Server

# Set pandas options to display all columns
pd.set_option('display.max_columns', None)


def main() -> None:
    """Main function to create a gene dictionary from Ensembl Biomart."""
    # Connect to Ensembl Biomart public server
    server = Server(host='http://www.ensembl.org')

    # Get the human dataset
    dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']

    # Query required fields
    df = dataset.query(attributes=[
        'external_gene_name',  # gene symbol
        'ensembl_gene_id',
        'ensembl_transcript_id',
        'chromosome_name',
        'start_position',
        'end_position',
        'transcript_mane_select'
    ])

    # Drop missing transcripts
    df = df.dropna(subset=['Transcript stable ID'])

    # Filter to standard chromosomes
    df = df[df['Chromosome/scaffold name'].isin([str(i) for i in range(1, 23)] + ['X', 'Y', 'MT'])]

    # Build dictionary
    final_dict = {}
    grouped = df.groupby('Gene name')

    for gene_symbol, group in grouped:
        gene_id = group['Gene stable ID'].iloc[0]
        all_transcripts = group['Transcript stable ID'].dropna().unique().tolist()

        mane_row = group.dropna(subset=['RefSeq match transcript (MANE Select)'])
        mane_transcript = mane_row['Transcript stable ID'].iloc[0] if not mane_row.empty else None

        # Genomic location from the first row (they're all from the same gene)
        loc = {
            'chromosome': group['Chromosome/scaffold name'].iloc[0],
            'start': int(group['Gene start (bp)'].min()),
            'end': int(group['Gene end (bp)'].max())
        }

        # Ensure mane_transcript is in the list
        if mane_transcript and mane_transcript not in all_transcripts:
            all_transcripts.append(mane_transcript)

        final_dict[gene_symbol] = {
            'gene_name': gene_id,
            'mane_transcript': mane_transcript,
            'all_transcripts': all_transcripts,
            'genomic_location': loc
        }

    # Export to JSON
    with open("final_dict_public.json", "w") as f:
        json.dump(final_dict, f, indent=2)

    print("Gene dictionary has been successfully created and saved to 'final_dict_public.json'.")


if __name__ == "__main__":
    main()
