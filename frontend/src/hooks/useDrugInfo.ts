import { useEffect, useState } from "react";
import axios from "axios";

import { backendUrl, searchTypes } from "../constants";
import { DrugProperties, URLParams } from "../types";

const toSubscript = (str: string): string => {
  const subscriptMap: { [key: string]: string } = {
    "0": "₀",
    "1": "₁",
    "2": "₂",
    "3": "₃",
    "4": "₄",
    "5": "₅",
    "6": "₆",
    "7": "₇",
    "8": "₈",
    "9": "₉",
  };
  return str
    .split("")
    .map((char) => subscriptMap[char] || char)
    .join("");
};

const useDrugInfo = (params: URLParams) => {
  const [drugInfo, setDrugInfo] = useState<DrugProperties[]>([]);
  const [error, setError] = useState<unknown | boolean>(false);


  const fetchDrugMolecule = async (drug: string) => {
    try {
      const response = await axios.get(
        `${backendUrl}/drug/get_molecule?drug_name=${drug}`,
        { responseType: "arraybuffer" },
      );
      const blob = new Blob([response.data], { type: "image/png" });
      return URL.createObjectURL(blob);
    } catch (error) {
      console.warn("Error fetching drug molecule:", error);
      return;
    }
  };

  useEffect(() => {
    setError(false);

    if (params.searchBy === searchTypes[2]?.param) {
      // When searching by side effect, we don't need to fetch drug info
      // Displaying only the side effect name
      setDrugInfo(
        params.terms.map((name) => ({
          name: name.charAt(0).toUpperCase() + name.slice(1),
        })),
      );
    } else {
      const getDrugInfo = async () => {
        try {
          const promises = params.terms.map((name) =>
            axios.get(
              `${backendUrl}/drug/get_info?drug_name=${name}&search_type=${params.searchBy}`,
            ),
          );
          const responses = await Promise.all(promises);

          const infoPromises = responses.flatMap((response) => {
            if (response.data.length > 0) {
              return response.data.map(async (drug: DrugProperties) => ({
                name: drug.name || "",
                classification: drug.classification || "",
                groups: drug.groups || "",
                iupac: drug.iupac || "",
                formula: drug.formula ? toSubscript(drug.formula) : "",
                brands: drug.brands || "",
                half_life: drug.half_life || "",
                indication: drug.indication || "",
                product: drug.product || "",
                molecule_url: (await fetchDrugMolecule(drug.name)) || "",
              })) as DrugProperties;
            } else {
              setError(true);
              return [];
            }
          });
          const infoArray: DrugProperties[] = await Promise.all(infoPromises);
          setDrugInfo(infoArray);
        } catch (error) {
          setError(true);
          setDrugInfo([]);
        }
      };

      void getDrugInfo();
    }
  }, [params.searchBy, params.terms]);

  return { drugInfo, drugInfoError: error };
};

export default useDrugInfo;
