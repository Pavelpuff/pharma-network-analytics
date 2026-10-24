import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import plotly.express as px


class PharmaVisualizer:
    """
    Advanced visualization module for pharmacological networks and EDA.
    """

    @staticmethod
    def plot_synonym_network(drugbank_id: str, df_synonyms: pd.DataFrame) -> None:
        """Plots a star-graph network of a drug and its associated synonyms."""
        G = nx.Graph()
        drug_data = df_synonyms[df_synonyms["DrugBank_ID"] == drugbank_id]
        if drug_data.empty or not drug_data["Synonyms"].values[0]:
            print(f"No synonyms found for DrugBank ID: {drugbank_id}")
            return

        drug_name = drug_data["Name"].values[0]
        G.add_node(drug_name, color="#e74c3c", size=1500)

        for synonym in drug_data["Synonyms"].values[0]:
            G.add_node(synonym, color="#3498db", size=800)
            G.add_edge(drug_name, synonym)

        pos = nx.spring_layout(G, seed=42)
        colors = [G.nodes[n]["color"] for n in G.nodes]
        sizes = [G.nodes[n]["size"] for n in G.nodes]

        fig, ax = plt.subplots(figsize=(12, 8), facecolor="white")
        ax.set_facecolor("white")

        nx.draw_networkx_nodes(
            G, pos, node_size=sizes, node_color=colors, alpha=0.85, ax=ax
        )
        nx.draw_networkx_edges(G, pos, width=1.5, edge_color="#bdc3c7", ax=ax)

        bbox_props = dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.9)
        nx.draw_networkx_labels(
            G, pos, font_size=10, font_weight="bold", bbox=bbox_props, ax=ax
        )

        plt.margins(0.2)
        plt.title(
            f"Synonym Network for {drug_name} (ID: {drugbank_id})",
            fontsize=16,
            pad=15,
            color="black",
        )
        plt.axis("off")
        plt.tight_layout()
        plt.show()

    @staticmethod
    def plot_bipartite_pathway_network(df_pathways: pd.DataFrame) -> None:
        """Plots a bipartite graph showing interactions between pathways and drugs."""
        G = nx.Graph()
        pathways = df_pathways["Pathway_Name"].unique()
        drugs = df_pathways["Drug_Name"].unique()

        G.add_nodes_from(pathways, bipartite=0, color="#e74c3c")
        G.add_nodes_from(drugs, bipartite=1, color="#3498db")

        edges = [
            (row["Pathway_Name"], row["Drug_Name"]) for _, row in df_pathways.iterrows()
        ]
        G.add_edges_from(edges)

        pos = nx.bipartite_layout(G, pathways)
        colors = [G.nodes[n]["color"] for n in G.nodes]

        fig, ax = plt.subplots(figsize=(14, 10), facecolor="white")
        ax.set_facecolor("white")

        nx.draw(
            G,
            pos,
            with_labels=True,
            node_color=colors,
            node_size=600,
            font_size=8,
            font_weight="bold",
            edge_color="#95a5a6",
            alpha=0.8,
            ax=ax,
        )

        plt.margins(0.15)
        plt.title(
            "Bipartite Network: Drug-Pathway Interactions", fontsize=16, color="black"
        )
        plt.axis("off")
        plt.show()

    @staticmethod
    def plot_pathway_counts(df_pathways: pd.DataFrame) -> None:
        """Generates a histogram of pathways per interacting drug."""
        pathway_counts = df_pathways["Drug_Name"].value_counts().reset_index()
        pathway_counts.columns = ["Drug_Name", "Pathway_Count"]
        pathway_counts = pathway_counts.sort_values(by="Pathway_Count", ascending=False)

        fig = px.bar(
            pathway_counts,
            x="Drug_Name",
            y="Pathway_Count",
            title="Total Pathways Interacting per Drug",
            labels={
                "Drug_Name": "Active Ingredient",
                "Pathway_Count": "Number of Pathways",
            },
            color_discrete_sequence=["skyblue"],
        )
        fig.update_layout(xaxis_tickangle=-45)
        fig.show()

    @staticmethod
    def plot_target_cellular_distribution(df_targets: pd.DataFrame) -> None:
        """Generates a pie chart of target cellular locations."""
        location_counts = df_targets["Cellular_Location"].value_counts().reset_index()
        location_counts.columns = ["Cellular_Location", "Count"]

        fig = px.pie(
            location_counts,
            names="Cellular_Location",
            values="Count",
            title="Target Distribution by Cellular Localization",
            color_discrete_sequence=px.colors.qualitative.Set3,
        )
        fig.show()

    @staticmethod
    def plot_drug_approval_status(df_statuses: pd.DataFrame) -> None:
        """Generates a pie chart of drug approval groups."""
        status_counts = df_statuses["Status_Group"].value_counts().reset_index()
        status_counts.columns = ["Status_Group", "Count"]

        fig = px.pie(
            status_counts,
            names="Status_Group",
            values="Count",
            title="Distribution of Pharmaceutical Approval Statuses",
            color_discrete_sequence=px.colors.qualitative.Set3,
        )
        fig.show()

    @staticmethod
    def plot_drug_interactions_count(df_interactions: pd.DataFrame) -> None:
        """Generates a histogram of drug-drug interactions per drug."""
        interaction_counts = df_interactions["Drug_A_Name"].value_counts().reset_index()
        interaction_counts.columns = ["Drug_Name", "Interaction_Count"]

        fig = px.bar(
            interaction_counts,
            x="Drug_Name",
            y="Interaction_Count",
            title="Number of Drug-Drug Interactions (DDI) per Drug",
            labels={
                "Drug_Name": "Drug Entity",
                "Interaction_Count": "Total Interactions",
            },
            color_discrete_sequence=["#636EFA"],
        )
        fig.update_layout(xaxis_tickangle=-45)
        fig.show()

    @staticmethod
    def plot_gene_drug_product_network(
        gene_name: str, df_targets: pd.DataFrame, df_products: pd.DataFrame
    ) -> None:
        """
        Custom hierarchical network mapping a specific Gene to interacting Drugs,
        and linking those Drugs to market Products.
        """
        gene_data = df_targets[df_targets["Gene_Name"] == gene_name]
        if gene_data.empty:
            print(f"No target data available for Gene: {gene_name}")
            return

        drug_ids = gene_data["DrugBank_ID"].unique()
        drugs_data = df_products[df_products["DrugBank_ID"].isin(drug_ids)]

        G = nx.Graph()
        G.add_node(gene_name, color="#e74c3c", size=1500)

        for _, row in gene_data.iterrows():
            drug_id = row["DrugBank_ID"]
            G.add_node(drug_id, color="#3498db", size=1000)
            G.add_edge(gene_name, drug_id)

        for _, row in drugs_data.iterrows():
            product_name = row["Product_Name"]
            drug_id = row["DrugBank_ID"]
            G.add_node(product_name, color="#2ecc71", size=600)
            G.add_edge(drug_id, product_name)

        fig, ax = plt.subplots(figsize=(16, 12), facecolor="white")
        ax.set_facecolor("white")

        colors = [G.nodes[n].get("color", "gray") for n in G.nodes]
        sizes = [G.nodes[n].get("size", 300) for n in G.nodes]

        pos = nx.spring_layout(G, k=0.9, seed=42)

        nx.draw_networkx_nodes(
            G, pos, node_color=colors, node_size=sizes, alpha=0.9, ax=ax
        )
        nx.draw_networkx_edges(G, pos, edge_color="#bdc3c7", alpha=0.7, ax=ax)

        bbox_props = dict(boxstyle="round,pad=0.2", fc="white", ec="gray", alpha=0.85)
        nx.draw_networkx_labels(
            G, pos, font_size=9, font_weight="bold", bbox=bbox_props, ax=ax
        )

        plt.margins(0.2)
        plt.title(
            f"Hierarchical Target Network: Gene ({gene_name}) → Active Drugs → Market Products",
            fontsize=16,
            pad=15,
            color="black",
        )
        plt.axis("off")
        plt.tight_layout()
        plt.show()
