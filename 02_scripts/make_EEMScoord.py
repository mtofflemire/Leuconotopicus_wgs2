import pandas as pd

infile = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/01_data/Leuconotopicus_albolarvatus_Meta.csv"
outfile = "/Users/michaeltofflemire/Mtofflemire Dropbox/Michael Tofflemire/Projects/Leuconotopicus_wgs/02_scripts/EEMS/eems.Leuconotopicus.coord"

df = pd.read_csv(infile)

df[["Longitude", "Latitude"]].to_csv(
    outfile,
    sep="\t",
    index=False,
    header=False
)