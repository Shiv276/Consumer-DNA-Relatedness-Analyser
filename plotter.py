import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

def predict_relationship(row):

    # Relationship criteria is based on the original KING paper
    # Classification boundaries are the geometric means between expected kinship coefficients, following conventional KING relationship thresholds
    #Manichaikul, A., Mychaleckyj, J. C., Rich, S. S., Daly, K., Sale, M., & Chen, W. M. (2010). Robust relationship inference in genome-wide association studies. Bioinformatics (Oxford, England), 26(22), 2867–2873. https://doi.org/10.1093/bioinformatics/btq559
    DUPLICATE_MIN = 1 / (2 ** (3 / 2))
    FIRST_DEGREE_MIN = 1 / (2 ** (5 / 2))
    SECOND_DEGREE_MIN = 1 / (2 ** (7 / 2))
    THIRD_DEGREE_MIN = 1 / (2 ** (9 / 2))

    #IBS0 val for Parent-Child relation is also published in a separate journal. Ideally IBS0 would be 0 in the absence of genotyping errors, however that is experimentally unviable.
    #0.0012 comes from observed IBS information in UK Biobank SNP-array data
    #Hofmeister RJ, Rubinacci S, Ribeiro DM, Buil A, Kutalik Z, Delaneau O. Parent-of-Origin inference for biobanks. Nat Commun. 2022 Nov 5;13(1):6668. doi: 10.1038/s41467-022-34383-6. PMID: 36335127; PMCID: PMC9637181.
    PARENT_CHILD_IBS0_MAX = 0.0012

    kinship = row['KINSHIP']
    ibs0 = row['IBS0']

    if kinship >= DUPLICATE_MIN:
        return "Duplicate / identical twin"

    if kinship >= FIRST_DEGREE_MIN:
        if ibs0 < PARENT_CHILD_IBS0_MAX:
            return "Parent-child-like"
        return "Full-sibling-like"

    if kinship >= SECOND_DEGREE_MIN:
        return "Second-degree"

    if kinship >= THIRD_DEGREE_MIN:
        return "Third-degree"

    return "Unrelated / distant"




def plotter(kin_file):
    df = pd.read_table(kin_file, sep=r"\s+") 

    # Form pair names and calculate predicted relation
    df["Pair"] = df["#FID1"].astype(str) + " vs " + df["FID2"].astype(str)
    df['Relation'] = df.apply(predict_relationship, axis=1)


    def make_label(row):
        # Show IBS0 for first-degree relationships because it is what separates parent-child-like from full-sibling-like
        if row["Relation"] in ["Parent-child-like", "Full-sibling-like"]:
            return (
                f"{row['Relation']} "
                f"(KIN={row['KINSHIP']:.3f}, IBS0={row['IBS0']:.4f})"
            )

        return f"{row['Relation']} (KIN={row['KINSHIP']:.3f})"

    df["Plot_label"] = df.apply(make_label, axis=1)


    #Show stronger relations at the top of df
    df = df.sort_values(by="KINSHIP", ascending=False)





    #Different colours for each relation
    colour_map = {
            "Duplicate / identical twin": "#9467bd",
            "Parent-child-like": "#1f77b4",
            "Full-sibling-like": "#ff7f0e",
            "Second-degree": "#2ca02c",
            "Third-degree": "#bcbd22",
            "Unrelated / distant": "#777777"
        }

    colours = df["Relation"].map(colour_map)


    fig, ax = plt.subplots(figsize=(14, 7))

    bars = ax.barh(
    df["Pair"],
    df["KINSHIP"],
    color=colours,
    edgecolor="black",
    height=0.8
    )

    ax.invert_yaxis() #Plot places first row at the bottom. Inversion allows stronger relations to be at the top while retaining proper DF order.

    # KING degree boundaries
    cutoffs = {
    "Third-degree cutoff": (0.0442, "#919191"),
    "Second-degree cutoff": (0.0884, "#424242"),
    "First-degree cutoff": (0.177, "#000000")
    }

    #Adds vertical line at degree boundaries based on cutoffs
    for name, (cutoff, colour) in cutoffs.items():
        ax.axvline(
            cutoff,
            color=colour,
            linestyle="--",
            linewidth=1.2,
            alpha=0.8,
            label=name
        )

    ax.legend(loc="lower right", frameon=False)


    # Add prediction labels beside each bar
    for y, row in enumerate(df.itertuples(index=False)): #itertuples() gives each DataFrame row, allowing columns to be accessed. y is index and row is each df row.

        text_x = max(row.KINSHIP, 0) + 0.003 #Moves label slightly to the right of the edge of each bar.
        ax.text(
            text_x,
            y,
            row.Plot_label,
            va="center",
            fontsize=10
        )

    # Leave enough space for negative bars and right-side labels
    left_limit = min(-0.03, df["KINSHIP"].min() - 0.01) #If the lowest KINSHIP goes below 0, set the left axis limit slightly below that
    right_limit = max(0.32, df["KINSHIP"].max() + 0.11) #Set the right axis limit higher than the max KINSHIP value to allow enough space for the longer relationship labels
    #These numbers were chosen arbitrarily to allow enough visual space between labels and axis limits
    ax.set_xlim(left_limit, right_limit)

    #Title for the plot
    ax.set_title(
        "Pairwise relatedness inferred from KING\n",
        fontsize=17
    )

    #Axis Ttiles. Y-title is not necessary.
    ax.set_xlabel("KING kinship coefficient", fontsize=12)
    ax.set_ylabel("")

    #One single line at the 0 mark rather than gridlines for a cleaner look.
    ax.axvline(
        x=0,
        color="black",
        linestyle="-",
        linewidth=1.2
    )

    plt.tight_layout()

    #Making a user-friendly version of the df in CSV format 
    user_df = df[['Pair', 'Relation', 'KINSHIP', 'IBS0', 'NSNP']]

    return user_df, fig

