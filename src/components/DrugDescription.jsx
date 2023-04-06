import React, { useState, useEffect } from 'react';
import axios from 'axios';

const drugDescriptionUrl = 'https://api.fda.gov/drug/drugsfda.json?search=openfda.brand_name:'

const DrugDescription = ({ drugName }) => {
    const [drugInfo, setDrugInfo] = useState({});

    useEffect(() => {
        const fetchDrugInfo = async () => {
            try {
                const response = await axios.get(
                    `${drugDescriptionUrl}${drugName}+openfda.generic_name:${drugName}&limit=1`
                );

                if (response.data.results && response.data.results.length > 0) {
                    const result = response.data.results[0];
                    console.log(result)
                    const info = {
                        drug_class: result.openfda.pharm_class_epc[0],
                        route: result.openfda.route[0],
                    };
                    setDrugInfo(info);
                } else {
                    setDrugInfo({ error: 'No information found.' });
                }
            } catch (error) {
                console.error('Error fetching drug information:', error);
                setDrugInfo({ error: 'Error fetching drug information.' });
            }
        };
        if (drugName) {
            fetchDrugInfo();
        }
    }, [drugName]);

    return (
        // align text to the left
        <div>
            {drugInfo.error ? (
                <p>{drugInfo.error}</p>
            ) : (
                <>
                    <p>Route: {drugInfo.route}</p>
                    <p>Drug Class: {drugInfo.drug_class}</p>
                </>
            )}
        </div>
    );
};

export default DrugDescription;
