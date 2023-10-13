from Bio import Entrez
from collections import defaultdict

# Provide your email address to NCBI
Entrez.email = "tbrunoir@ucdavis.edu"

# Read genera from tmpA.txt
with open("tmpA", "r") as file:
    genera = [line.strip() for line in file]

# Create a dictionary to store unique genera and corresponding taxids
genus_taxid_dict = defaultdict(set)

# Fetch TaxIDs for each genus
for genus in genera:
    try:
        handle = Entrez.esearch(db="taxonomy", term=f"{genus} [ORGN]", retmode="json")
        record = Entrez.read(handle)
        taxids = record["esearchresult"]["idlist"]
        for taxid in taxids:
            genus_taxid_dict[genus].add(taxid)
    except Exception as e:
        print(f"Error fetching TaxIDs for {genus}: {e}")

# Write the genus and taxid information to tmpB_taxids
with open("tmpB_taxids", "w") as output_file:
    output_file.write("Genus,TaxID\n")  # Write header
    for genus, taxids in genus_taxid_dict.items():
        for taxid in taxids:
            output_file.write(f"{genus},{taxid}\n")
