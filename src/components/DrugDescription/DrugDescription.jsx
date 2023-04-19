import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { CapsulePill as CapsulePillIcon } from "react-bootstrap-icons";
import { Card, Placeholder } from 'react-bootstrap';
import './DrugDescription.css';

const drugDescriptionUrl = 'https://api.fda.gov/drug/drugsfda.json?search=openfda.brand_name:'

const DrugDescription = ({ drugName }) => {
    const [drugInfo, setDrugInfo] = useState({});
    const [isLoading, setIsLoading] = useState(true);

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
                        drug_name: result.openfda.substance_name[0].charAt(0).toUpperCase() + result.openfda.substance_name[0].slice(1).toLowerCase(),
                        drug_class: result.openfda.pharm_class_epc[0].replace(/\[EPC]/, ''),
                        route: result.openfda.route[0].charAt(0).toUpperCase() + result.openfda.route[0].slice(1).toLowerCase(),
                    };
                    setDrugInfo(info);
                } else {
                    setDrugInfo({ error: 'No information found.' });
                }
            } catch (error) {
                console.error('Error fetching drug information:', error);
                setDrugInfo({ error: 'No information found.' });
            }
        };
        if (drugName) {
            setIsLoading(true);
            fetchDrugInfo();
            setIsLoading(false);
        }
    }, [drugName]);

    return (
        <Card>
            <Card.Header><CapsulePillIcon className="drug-badge mx-1" />
                {isLoading ? (
                    <Placeholder as="span" animation="glow">
                        <Placeholder xs={4} />
                    </Placeholder>
                ) : (
                    drugInfo.drug_name
                )}
            </Card.Header>
            {drugInfo.error ? (
                <p>{drugInfo.error}</p>
            ) : (
                <Card.Body className="text-md-start">
                    <p>Administration:&nbsp;
                        {isLoading ? (
                            <Placeholder as="span" animation="glow">
                                <Placeholder xs={4} />
                            </Placeholder>
                        ) : (
                            drugInfo.route
                        )}
                    </p>
                    <p>Drug Class:&nbsp;
                        {isLoading ? (
                            <Placeholder as="span" animation="glow">
                                <Placeholder xs={4} />
                            </Placeholder>
                        ) : (
                            drugInfo.drug_class
                        )}
                    </p>
                </Card.Body>
            )}
        </Card>
    );
};

export default DrugDescription;
