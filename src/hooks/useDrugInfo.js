import {useEffect, useState} from "react";
import axios from "axios";

function toSubscript(str) {
    return str.replace(/\d/g, digit => '₀₁₂₃₄₅₆₇₈₉'[digit]);
}

const useDrugInfo = (drugName, searchType) => {
    const drugAPI = 'http://localhost:16000/api'

    const [drugInfo, setDrugInfo] = useState([])
    const [isLoading, setIsLoading] = useState(false)

    useEffect(() => {
        const getDrugInfo = async () => {
            setIsLoading(true)
            try {
                const promises = drugName.map(name =>
                    axios.get(`${drugAPI}/get_info?drug_name=${name}&search_type=${searchType}`)
                );
                const responses = await Promise.all(promises);

                const infoPromises = responses.flatMap(response => {
                    if (response.data.length > 0) {
                        return response.data.map(async drug => ({
                            drug_name: drug.name,
                            drug_class: drug.classification,
                            groups: drug.groups,
                            iupac: drug.iupac,
                            formula: toSubscript(drug.formula),
                            brands: drug.brands,
                            half_life: drug.half_life,
                            indication: drug.indication,
                            product: drug.product,
                            molecule_url: await fetchDrugMolecule(drug.name)
                        }));
                    } else {
                        return [];
                    }
                });
                const infoArray = await Promise.all(infoPromises);
                setDrugInfo(infoArray);
            } catch (error) {
                console.error('Error fetching drug information:', error);
            }
        };

        const fetchDrugMolecule = async (drugInfo) => {
            try {
                const response = await axios.get(
                    `${drugAPI}/get_molecule?drug_name=${drugInfo}`,
                    { responseType: 'arraybuffer' }
                )
                const blob = new Blob([response.data], { type: 'image/png' })
                return URL.createObjectURL(blob)
            } catch (error) {
                console.error('Error fetching drug molecule:', error)
            }
        }

        getDrugInfo()
    }, [drugName, searchType])

    return { drugInfo, isLoading }
}

export default useDrugInfo