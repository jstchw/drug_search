import {useEffect, useState} from "react";
import axios from "axios";
import {searchTypes} from "../components/OptionModal/OptionModal";

const toSubscript = (str) => {
    return str.replace(/\d/g, digit => '₀₁₂₃₄₅₆₇₈₉'[digit]);
}

const useDrugInfo = (drugName, searchType) => {
    const drugAPI = window.REACT_APP_API_URL

    const [drugInfo, setDrugInfo] = useState([])

    const fetchDrugMolecule = async (drugInfo) => {
        try {
            const response = await axios.get(
                `${drugAPI}/get_molecule?drug_name=${drugInfo}`,
                { responseType: 'arraybuffer' }
            )
            const blob = new Blob([response.data], { type: 'image/png' })
            return URL.createObjectURL(blob)
        } catch (error) {
            console.warn('Error fetching drug molecule:', error)
        }
    }

    useEffect(() => {
        if(searchType === searchTypes[2].value) {
            setDrugInfo(drugName.map(name => ({
                ADE: name.charAt(0).toUpperCase() + name.slice(1),
            })))
        } else {
            const getDrugInfo = async () => {
                try {
                    const promises = drugName.map(name =>
                        axios.get(`${drugAPI}/get_info?drug_name=${name}&search_type=${searchType}`)
                    );
                    const responses = await Promise.all(promises);

                    const infoPromises = responses.flatMap(response => {
                        if (response.data.length > 0) {
                            return response.data.map(async drug => ({
                                drug_name: drug.name || "",
                                drug_class: drug.classification || "",
                                groups: drug.groups || "",
                                iupac: drug.iupac || "",
                                formula: drug.formula ? toSubscript(drug.formula) : "",
                                brands: drug.brands || "",
                                half_life: drug.half_life || "",
                                indication: drug.indication || "",
                                product: drug.product || "",
                                molecule_url: await fetchDrugMolecule(drug.name) || "",
                                ADE: null,
                                full_info: true,
                            }));
                        } else {
                            return drugName.map(name => ({
                                drug_name: name.charAt(0).toUpperCase() + name.slice(1),
                                ADE: null,
                                full_info: false,
                            }));
                        }
                    });
                    const infoArray = await Promise.all(infoPromises);
                    setDrugInfo(infoArray);
                } catch (error) {
                    console.warn('Error fetching drug information:', error);
                }
            };

            getDrugInfo()
        }

    }, [JSON.stringify(drugName)]) // eslint-disable-line react-hooks/exhaustive-deps

    return drugInfo
}

export default useDrugInfo