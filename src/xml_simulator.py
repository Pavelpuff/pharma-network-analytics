import random
import xml.etree.ElementTree as ET


def generate_mock_database(
    input_path: str, output_path: str, total_records: int = 20000
) -> None:
    """
    Parses the original XML, extracts attribute pools from existing drugs,
    and generates synthetic drug records to reach the total_records target.
    """
    ET.register_namespace("", "http://www.drugbank.ca")
    tree = ET.parse(input_path)
    root = tree.getroot()
    ns = {"db": "http://www.drugbank.ca"}

    existing_drugs = root.findall("db:drug", namespaces=ns)
    if not existing_drugs:
        raise ValueError("No drugs found in the source XML.")

    tag_pools = {}
    tags_to_extract = [
        "name",
        "description",
        "indication",
        "pharmacodynamics",
        "mechanism-of-action",
        "toxicity",
    ]

    for tag in tags_to_extract:
        tag_pools[tag] = [
            drug.findtext(f"db:{tag}", default="", namespaces=ns)
            for drug in existing_drugs
            if drug.find(f"db:{tag}", namespaces=ns) is not None
        ]

    max_id_num = 0
    for drug in existing_drugs:
        primary_id = drug.findtext(
            "db:drugbank-id[@primary='true']", default="", namespaces=ns
        )
        if primary_id.startswith("DB"):
            try:
                num = int(primary_id.replace("DB", ""))
                max_id_num = max(max_id_num, num)
            except ValueError:
                continue

    drugs_to_generate = total_records - len(existing_drugs)

    print(f"Generating {drugs_to_generate} synthetic drug records...")

    for i in range(1, drugs_to_generate + 1):
        new_drug = ET.Element("{http://www.drugbank.ca}drug")
        new_drug.set("type", random.choice(["biotech", "small molecule"]))

        new_id_val = f"DB{str(max_id_num + i).zfill(5)}"
        id_elem = ET.SubElement(
            new_drug, "{http://www.drugbank.ca}drugbank-id", primary="true"
        )
        id_elem.text = new_id_val

        for tag, values in tag_pools.items():
            if values:
                elem = ET.SubElement(new_drug, f"{{http://www.drugbank.ca}}{tag}")
                elem.text = random.choice(values)

        root.append(new_drug)

    tree.write(output_path, encoding="utf-8", xml_declaration=True)
    print(
        f"Mock database saved to {output_path} with {len(root.findall('db:drug', namespaces=ns))} records."
    )


if __name__ == "__main__":
    generate_mock_database(
        "data/drugbank_partial.xml", "data/drugbank_partial_and_generated.xml"
    )
