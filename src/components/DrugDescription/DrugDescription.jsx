import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { CapsulePill as CapsulePillIcon } from "react-bootstrap-icons";
import { Card, Placeholder } from 'react-bootstrap';
import './DrugDescription.css';

const drugDescriptionUrl = 'http://localhost:16000/api';

const DrugDescription = ({ drugName }) => {

    const [drugInfo, setDrugInfo] = useState({});
    const [isLoading, setIsLoading] = useState(true);

    const [moleculeUrl, setMoleculeUrl] = useState(null);

    const fetchDrugInfo = async () => {
        try {
            const response = await axios.get(
                `${drugDescriptionUrl}/get_info?drug_name=${drugName}`
            );

            // Data from API is managed here
            if (response.data.length > 0) {
                const info = {
                    drug_name: response.data[0].name,
                    drug_class: response.data[0].classification
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

    const fetchDrugMolecule = async () => {
        try {
            const response = await axios.get(
                `${drugDescriptionUrl}/get_molecule?drug_name=${drugName}`,
                { responseType: 'arraybuffer' }
            )
            const blob = new Blob([response.data], { type: 'image/png' })
            const url = URL.createObjectURL(blob)
            setMoleculeUrl(url)
        } catch (error) {
            console.error('Error fetching drug molecule:', error)
        }
    }

    useEffect(() => {
        if (drugName) {
            setIsLoading(true);
            fetchDrugMolecule();
            fetchDrugInfo();
            setIsLoading(false);
        }
    }, [drugName]);

    return (
        <Card>
            <Card.Header><CapsulePillIcon className="drug-badge mx-1" />
                <span>General Knowledge</span>
            </Card.Header>
            <Card.Img src={moleculeUrl} alt={drugName} className="molecule" style={{ backgroundColor: 'transparent'}} />

            {drugInfo.error ? (
                <p>{drugInfo.error}</p>
            ) : (
                <Card.Body className="text-md-start">
                    {/*<p>Half-Life:&nbsp;*/}
                    {/*    {isLoading ? (*/}
                    {/*        <Placeholder as="span" animation="glow">*/}
                    {/*            <Placeholder xs={4} />*/}
                    {/*        </Placeholder>*/}
                    {/*    ) : (*/}
                    {/*        drugInfo.half_life*/}
                    {/*    )}*/}
                    {/*</p>*/}
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
