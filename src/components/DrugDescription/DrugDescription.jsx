import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { CapsulePill as CapsulePillIcon } from "react-bootstrap-icons";
import { Card, Placeholder } from 'react-bootstrap';
import './DrugDescription.css';

const drugDescriptionUrl = 'http://localhost:16000/api/get_info';

const DrugDescription = ({ drugName }) => {
    const [drugInfo, setDrugInfo] = useState({});
    const [isLoading, setIsLoading] = useState(true);

    useEffect(() => {
        const fetchDrugInfo = async () => {
            try {
                const response = await axios.get(
                    `${drugDescriptionUrl}?drug_name=${drugName}`
                );

                if (response.data.length > 0) {
                    const info = {
                        drug_name: response.data[0].name,
                        drug_class: response.data[0].classification,
                        half_life: response.data[0].half_life,
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
                    <p>Half-Life:&nbsp;
                        {isLoading ? (
                            <Placeholder as="span" animation="glow">
                                <Placeholder xs={4} />
                            </Placeholder>
                        ) : (
                            drugInfo.half_life
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
