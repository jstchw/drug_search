import React from 'react';
import {Capsule as CapsuleIcon} from "react-bootstrap-icons";
import { Card, Placeholder, Badge } from 'react-bootstrap';
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
