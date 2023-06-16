import React, { useState, useEffect } from 'react';
import axios from 'axios';
import { Capsule as CapsuleIcon } from "react-bootstrap-icons";
import { Card, Placeholder, Badge } from 'react-bootstrap';
import './DrugDescription.css';

const drugDescriptionUrl = 'http://localhost:16000/api';

const DrugDescription = ({ drugName, onRetrieved }) => {
    const [drugInfo, setDrugInfo] = useState({});
    const [isLoading, setIsLoading] = useState(false);

    const [moleculeUrl, setMoleculeUrl] = useState(null);


    // Very weird function to convert the formula numbers to subscript
    function toSubscript(str) {
        return str.replace(/\d/g, digit => '₀₁₂₃₄₅₆₇₈₉'[digit]);
    }


    useEffect(() => {
        // Prevents the app from crashing if something is not available
        const handleRetrieved = (info) => {
            const data = Object.entries(info).reduce(
                (acc, [key, value]) => {
                    if (value !== null) {
                        acc[key] = value;
                    } else {
                        acc[key] = 'N/A';
                    }
                    return acc;
                }, {})
            onRetrieved(data);
        }

        const fetchDrugInfo = async () => {
            try {
                const response = await axios.get(
                    `${drugDescriptionUrl}/get_info?drug_name=${drugName}`
                );

                // Data from API is managed here
                if (response.data.length > 0) {
                    const info = {
                        drug_name: response.data[0].name,
                        drug_class: response.data[0].classification,
                        groups: response.data[0].groups,
                        iupac: response.data[0].iupac,
                        formula: toSubscript(response.data[0].formula),
                        brands: response.data[0].brands,
                        half_life: response.data[0].half_life,
                        indication: response.data[0].indication,
                    };
                    setDrugInfo(info);
                    // Pass information to parent component to display name and groups
                    handleRetrieved(info)
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
        const fetchAllData = async () => {
            if (drugName) {
                setIsLoading(true);
                await fetchDrugMolecule();
                await fetchDrugInfo();
                setIsLoading(false);
            }
        }

        fetchAllData()
    }, [drugName, onRetrieved]);

    return (
        <Card className={'drug-description-card'}>
            <Card.Header>
                <CapsuleIcon className="drug-badge mx-1" />
                <span>Molecule</span>
            </Card.Header>

            {moleculeUrl && <Card.Img src={moleculeUrl} alt={drugName} className="molecule" style={{ backgroundColor: 'transparent'}} />}

            <Card.Body className="text-md-start">
                <p><Badge>IUPAC</Badge>&nbsp;
                    {isLoading ? (
                        <>
                            <Placeholder as='span' animation="glow">
                                <Placeholder xs={4} />
                            </Placeholder>
                        </>
                    ) : (
                        drugInfo.iupac
                    )}
                </p>
                <p><Badge>Class</Badge>&nbsp;
                    {isLoading ? (
                        <>
                            <Placeholder as='span' animation="glow">
                                <Placeholder xs={4} />
                            </Placeholder>
                        </>
                    ) : (
                        drugInfo.drug_class
                    )}
                </p>
                <p><Badge>Formula</Badge>&nbsp;
                    {isLoading ? (
                        <>
                            <Placeholder as='span' animation="glow">
                                <Placeholder xs={4} />
                            </Placeholder>
                        </>
                    ) : (
                        drugInfo.formula
                    )}
                </p>
            </Card.Body>
        </Card>
    );
};

export default DrugDescription;
