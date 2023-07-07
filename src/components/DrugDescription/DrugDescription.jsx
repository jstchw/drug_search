import React from 'react';
import {Capsule as CapsuleIcon} from "react-bootstrap-icons";
import { Card, Badge } from 'react-bootstrap';
import './DrugDescription.css';


const DrugDescription = ({ drugInfo, showAdditionalSearch }) => {

    const dualSearchRef = React.useRef(showAdditionalSearch);

    return (
        <Card className={'drug-description-card mb-4'}>
            <Card.Header>
                <CapsuleIcon className="drug-badge mx-1" />
                {dualSearchRef.current ? (
                    <>
                        <span className={'drug-name'}>{drugInfo.drug_name}</span>
                    </>
                ) : (
                    <>
                        <span className={'drug-name'}>{'Molecule'}</span>
                    </>
                )}
            </Card.Header>

            {drugInfo.molecule_url && <Card.Img src={drugInfo.molecule_url} alt={drugInfo.drug_name}
                                                className="molecule" style={{ backgroundColor: 'transparent'}} />}

            <Card.Body className="text-md-start">
                {drugInfo && drugInfo.iupac &&
                    <p><Badge>IUPAC</Badge>&nbsp;<span className={'iupac'}>{drugInfo.iupac}</span></p>
                }
                {drugInfo && drugInfo.drug_class &&
                    <p><Badge>Class</Badge>&nbsp;{drugInfo.drug_class}</p>
                }
                {drugInfo && drugInfo.formula &&
                    <p><Badge>Formula</Badge>&nbsp;{drugInfo.formula}</p>
                }
            </Card.Body>
        </Card>
    );
};

export default DrugDescription;
