from selenium import webdriver
from selenium.webdriver.common.by import By
from lib_manejo_csv import escribe_csv, crea_csv, lee_csv


def get_uniprot_ids_from_complexes(file_path: str, file_path_out: str):
    _protens_list = lee_csv(file_path)
    crea_csv(file_path_out, ["db_id", "db_name", "db_uniprot_id"])
    _final_list = []
    _url_base = "https://www.yeastgenome.org/"
    _final_url = "/protein#external_ids"
    # Set up the webdriver
    driver = webdriver.Chrome()
    for _prot in _protens_list:
        try:
            _protein_object = {
                "db_id": _prot[0],
                "db_name": _prot[1],
                "db_uniprot_id": "",
            }
            driver.get(_url_base + "locus/" + _prot[1] + _final_url)
            externals = driver.find_element(By.ID, "external_ids")
            _search_input = externals.find_element(
                By.CLASS_NAME, "right"
            ).find_element(  # noqa: E501
                By.TAG_NAME, "input"  # noqa: E501
            )
            _search_input.send_keys("uniprot")
            _externals_id = externals.find_elements(By.TAG_NAME, "tr")
            for _id in _externals_id:
                _data = (_id.text).split(" ")
                if _data[1] == "UniProtKB":
                    _protein_object["db_uniprot_id"] = _data[0]
                    break
            _final_list.append(_protein_object)
        except Exception as e:
            escribe_csv("error_scrap_logs.csv", [_prot, e])
            pass
    # Close the webdriver
    driver.quit()
    # Write the results
    for _prot in _final_list:
        escribe_csv(
            file_path_out,
            [_prot["db_id"], _prot["db_name"], _prot["db_uniprot_id"]],
        )


# Main
if __name__ == "__main__":
    # for i in range(1, 17):
    i = 17
    _name = "./protein_db_chunk/proteins_all_db_chunk" + str(i) + ".csv"
    _out_name = (
        "./uniprot_mapping/proteins_all_db_chunk" + str(i) + "_uniprot.csv"
    )  # noqa: E501
    try:
        print("LOGS: Processing " + _name)
        get_uniprot_ids_from_complexes(_name, _out_name)
        print("LOGS: Finished " + _name)
    except Exception as e:
        escribe_csv("error_scrap_logs.csv", [_name, e])
        pass
