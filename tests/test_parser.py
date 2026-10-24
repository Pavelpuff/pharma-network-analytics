import pandas as pd
import pytest

from src.data_parser import DrugBankParser


@pytest.fixture(scope="module")
def setup_dummy_xml(tmp_path_factory):
    """Creates a temporary XML file for testing to avoid dependencies on real data."""
    xml_content = """<?xml version="1.0" encoding="UTF-8"?>
    <drugbank xmlns="http://www.drugbank.ca">
        <drug type="biotech">
            <drugbank-id primary="true">DB00001</drugbank-id>
            <name>Lepirudin</name>
            <description>Thrombin inhibitor</description>
            <pathways>
                <pathway>
                    <smpdb-id>SMP0000278</smpdb-id>
                    <name>Lepirudin Action Pathway</name>
                    <category>drug_action</category>
                    <drugs>
                        <drug>
                            <drugbank-id>DB00001</drugbank-id>
                            <name>Lepirudin</name>
                        </drug>
                        <drug>
                            <drugbank-id>DB01373</drugbank-id>
                            <name>Calcium</name>
                        </drug>
                    </drugs>
                </pathway>
            </pathways>
        </drug>
    </drugbank>
    """
    temp_dir = tmp_path_factory.mktemp("data")
    xml_file = temp_dir / "test_drugbank.xml"
    xml_file.write_text(xml_content, encoding="utf-8")
    return str(xml_file)


def test_extract_basic_info(setup_dummy_xml):
    parser = DrugBankParser(setup_dummy_xml)
    df = parser.extract_basic_info()

    assert isinstance(df, pd.DataFrame)
    assert not df.empty
    assert df.iloc[0]["DrugBank_ID"] == "DB00001"
    assert df.iloc[0]["Name"] == "Lepirudin"


def test_extract_pathways(setup_dummy_xml):
    parser = DrugBankParser(setup_dummy_xml)
    df = parser.extract_pathways()

    assert isinstance(df, pd.DataFrame)
    assert len(df) == 2  # Lepirudin + Calcium in the pathway
    assert "Calcium" in df["Drug_Name"].values
    assert df.iloc[0]["SMPDB_ID"] == "SMP0000278"
