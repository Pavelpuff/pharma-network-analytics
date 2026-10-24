import pandas as pd
from lxml import etree


class DrugBankParser:
    """
    Memory-efficient XML parser for the DrugBank database using lxml iterparse.
    Extracts nested pharmacological data into structured pandas DataFrames.
    """

    def __init__(self, file_path: str):
        self.file_path = file_path
        self.ns = {"db": "http://www.drugbank.ca"}

    def _get_text(self, element: etree.Element, xpath: str, default: str = "") -> str:
        found = element.find(xpath, namespaces=self.ns)
        return found.text if found is not None and found.text else default

    def _clear_element(self, elem: etree.Element) -> None:
        """Clears the element to free memory during iterative parsing."""
        elem.clear()
        while elem.getprevious() is not None:
            del elem.getparent()[0]

    def extract_basic_info(self) -> pd.DataFrame:
        """Extracts basic drug info: ID, Name, Type, Description, Form, Indication, Mechanism, Food Interactions."""
        data = []
        for _, elem in etree.iterparse(
            self.file_path, events=("end",), tag=f"{{{self.ns['db']}}}drug"
        ):
            drug_id = self._get_text(elem, "db:drugbank-id[@primary='true']")
            if drug_id:
                food_interactions = [
                    fi.text
                    for fi in elem.findall(
                        "db:food-interactions/db:food-interaction", namespaces=self.ns
                    )
                    if fi.text
                ]
                data.append(
                    {
                        "DrugBank_ID": drug_id,
                        "Name": self._get_text(elem, "db:name"),
                        "Type": elem.get("type", ""),
                        "Description": self._get_text(elem, "db:description"),
                        "Form": self._get_text(elem, "db:state"),
                        "Indication": self._get_text(elem, "db:indication"),
                        "Mechanism_of_Action": self._get_text(
                            elem, "db:mechanism-of-action"
                        ),
                        "Food_Interactions": ", ".join(food_interactions),
                    }
                )
            self._clear_element(elem)
        return pd.DataFrame(data)

    def extract_drug_synonyms(self) -> pd.DataFrame:
        """Extracts DrugBank IDs, Names, and Synonyms."""
        data = []
        for _, elem in etree.iterparse(
            self.file_path, events=("end",), tag=f"{{{self.ns['db']}}}drug"
        ):
            drug_id = self._get_text(elem, "db:drugbank-id[@primary='true']")
            if drug_id:
                synonyms = [
                    s.text
                    for s in elem.findall("db:synonyms/db:synonym", namespaces=self.ns)
                    if s.text
                ]
                data.append(
                    {
                        "DrugBank_ID": drug_id,
                        "Name": self._get_text(elem, "db:name"),
                        "Synonyms": synonyms,
                    }
                )
            self._clear_element(elem)
        return pd.DataFrame(data)

    def extract_products(self) -> pd.DataFrame:
        """Extracts pharmaceutical products associated with the drugs."""
        data = []
        for _, elem in etree.iterparse(
            self.file_path, events=("end",), tag=f"{{{self.ns['db']}}}drug"
        ):
            drug_id = self._get_text(elem, "db:drugbank-id[@primary='true']")
            if drug_id:
                for product in elem.findall(
                    "db:products/db:product", namespaces=self.ns
                ):
                    data.append(
                        {
                            "DrugBank_ID": drug_id,
                            "Product_Name": self._get_text(product, "db:name"),
                            "Manufacturer": self._get_text(product, "db:labeller"),
                            "NDC_Code": self._get_text(product, "db:ndc-product-code"),
                            "Dosage_Form": self._get_text(product, "db:dosage-form"),
                            "Route": self._get_text(product, "db:route"),
                            "Strength": self._get_text(product, "db:strength"),
                            "Country": self._get_text(product, "db:country"),
                            "Regulatory_Agency": self._get_text(product, "db:source"),
                        }
                    )
            self._clear_element(elem)
        return pd.DataFrame(data)

    def extract_pathways(self) -> pd.DataFrame:
        """Extracts pathways and all specific drugs interacting with them."""
        data = []
        for _, elem in etree.iterparse(
            self.file_path, events=("end",), tag=f"{{{self.ns['db']}}}drug"
        ):
            if elem.getparent().tag != f"{{{self.ns['db']}}}drugbank":
                continue

            drug_id = self._get_text(elem, "db:drugbank-id[@primary='true']")

            if drug_id:
                for pathway in elem.findall(
                    "db:pathways/db:pathway", namespaces=self.ns
                ):
                    pathway_name = self._get_text(pathway, "db:name")
                    category = self._get_text(pathway, "db:category")
                    smpdb_id = self._get_text(pathway, "db:smpdb-id")

                    for p_drug in pathway.findall(
                        "db:drugs/db:drug", namespaces=self.ns
                    ):
                        data.append(
                            {
                                "Pathway_Name": pathway_name,
                                "Category": category,
                                "SMPDB_ID": smpdb_id,
                                "Drug_Name": self._get_text(p_drug, "db:name"),
                                "DrugBank_ID": self._get_text(p_drug, "db:drugbank-id"),
                            }
                        )
            self._clear_element(elem)

        return pd.DataFrame(data).drop_duplicates(subset=["Pathway_Name", "Drug_Name"])

    def extract_drug_interactions(self) -> pd.DataFrame:
        """Extracts drug-drug interactions (DDI)."""
        data = []
        for _, elem in etree.iterparse(
            self.file_path, events=("end",), tag=f"{{{self.ns['db']}}}drug"
        ):
            drug_id = self._get_text(elem, "db:drugbank-id[@primary='true']")
            drug_name = self._get_text(elem, "db:name")
            if drug_id:
                for interaction in elem.findall(
                    "db:drug-interactions/db:drug-interaction", namespaces=self.ns
                ):
                    data.append(
                        {
                            "Drug_A_ID": drug_id,
                            "Drug_A_Name": drug_name,
                            "Drug_B_ID": self._get_text(interaction, "db:drugbank-id"),
                            "Drug_B_Name": self._get_text(interaction, "db:name"),
                        }
                    )
            self._clear_element(elem)
        return pd.DataFrame(data)

    def extract_targets(self) -> pd.DataFrame:
        """Extracts target proteins and their genetic/cellular attributes."""
        data = []
        for _, elem in etree.iterparse(
            self.file_path, events=("end",), tag=f"{{{self.ns['db']}}}drug"
        ):
            drug_id = self._get_text(elem, "db:drugbank-id[@primary='true']")
            if drug_id:
                for target in elem.findall("db:targets/db:target", namespaces=self.ns):
                    polypeptide = target.find("db:polypeptide", namespaces=self.ns)
                    if polypeptide is not None:
                        atlas_id, ext_id, source = "", "", ""
                        for ext in polypeptide.findall(
                            "db:external-identifiers/db:external-identifier",
                            namespaces=self.ns,
                        ):
                            resource = self._get_text(ext, "db:resource")
                            identifier = self._get_text(ext, "db:identifier")
                            if resource == "GenAtlas":
                                atlas_id = identifier
                            if not ext_id:
                                source, ext_id = resource, identifier

                        data.append(
                            {
                                "DrugBank_ID": drug_id,
                                "Target_ID": self._get_text(target, "db:id"),
                                "Target_Name": self._get_text(target, "db:name"),
                                "Source": source,
                                "External_ID": ext_id,
                                "Polypeptide_Name": self._get_text(
                                    polypeptide, "db:name"
                                ),
                                "Gene_Name": self._get_text(
                                    polypeptide, "db:gene-name"
                                ),
                                "GenAtlas_ID": atlas_id,
                                "Chromosome_Location": self._get_text(
                                    polypeptide, "db:chromosome-location"
                                ),
                                "Cellular_Location": self._get_text(
                                    polypeptide, "db:cellular-location"
                                ),
                            }
                        )
            self._clear_element(elem)
        return pd.DataFrame(data)

    def extract_drug_statuses(self) -> pd.DataFrame:
        """Extracts drug approval statuses."""
        data = []
        for _, elem in etree.iterparse(
            self.file_path, events=("end",), tag=f"{{{self.ns['db']}}}drug"
        ):
            drug_id = self._get_text(elem, "db:drugbank-id[@primary='true']")
            drug_name = self._get_text(elem, "db:name")
            if drug_id:
                groups = [
                    g.text
                    for g in elem.findall("db:groups/db:group", namespaces=self.ns)
                    if g.text
                ]
                for group in groups:
                    if group in ["experimental", "investigational"]:
                        group = "experimental/investigational"
                    data.append(
                        {
                            "DrugBank_ID": drug_id,
                            "Drug_Name": drug_name,
                            "Status_Group": group,
                        }
                    )
            self._clear_element(elem)
        return pd.DataFrame(data)
