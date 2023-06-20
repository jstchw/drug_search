import React, { useState, useEffect } from 'react';
import axios from 'axios';
import {Capsule as CapsuleIcon} from "react-bootstrap-icons";
import { Card, Placeholder, Badge } from 'react-bootstrap';
import './DrugDescription.css';

const drugDescriptionUrl = 'http://localhost:16000/api';

const DrugDescription = ({ drugName, onRetrieved, showAdditionalSearch }) => {
    const [drugInfo, setDrugInfo] = useState({});

    const [moleculeUrl, setMoleculeUrl] = useState(null);

    const prevDrugNameRef = React.useRef();

    const [dualSearch, setDualSearch] = useState(showAdditionalSearch);


    // Very weird function to convert the formula numbers to subscript
    function toSubscript(str) {
        return str.replace(/\d/g, digit => '₀₁₂₃₄₅₆₇₈₉'[digit]);
    }

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

        if(showAdditionalSearch) {
            setDualSearch(true)
            onRetrieved(prev => [...prev, data])
        } else {
            setDualSearch(false)
            onRetrieved([data])
        }
    }


    useEffect(() => {
        if(drugInfo.drug_name) {
            handleRetrieved(drugInfo)
        }
    }, [drugInfo])


    useEffect(() => {

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
                    }
                    setDrugInfo(info)
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
            if (drugName && drugName !== prevDrugNameRef.current) {
                try {
                    await Promise.all([fetchDrugMolecule(), fetchDrugInfo()])
                } catch (error) {
                    console.error('Error fetching all data:', error)
                }
            }
        }

        fetchAllData()
        prevDrugNameRef.current = drugName;
    }, [drugName])

    return (
        <Card className={'drug-description-card mb-4'}>
            <Card.Header>
                <CapsuleIcon className="drug-badge mx-1" />
                {dualSearch ? (
                    <span>{drugInfo.drug_name}</span>
                ) : (
                    <span>Molecule</span>
                )}
            </Card.Header>

            {moleculeUrl && <Card.Img src={moleculeUrl} alt={drugName} className="molecule" style={{ backgroundColor: 'transparent'}} />}

            <Card.Body className="text-md-start">
                <p><Badge>IUPAC</Badge>&nbsp;
                    {!drugInfo ? (
                        <>
                            <Placeholder as='span' animation="glow">
                                <Placeholder xs={4} />
                            </Placeholder>
                        </>
                    ) : (
                        <span className={'iupac'}>{drugInfo.iupac}</span>
                    )}
                </p>
                <p><Badge>Class</Badge>&nbsp;
                    {!drugInfo ? (
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
                    {!drugInfo ? (
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
