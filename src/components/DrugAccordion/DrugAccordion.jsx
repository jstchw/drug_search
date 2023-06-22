import {Accordion, Badge, Placeholder} from "react-bootstrap";
import React from "react";

const DrugAccordion = (props) => {
    return (
        <Accordion>
            <Accordion.Item eventKey={'0'}>
                <Accordion.Header><Badge>ADEs Reported</Badge></Accordion.Header>
                <Accordion.Body>
                    {parseInt(props.totalCount).toLocaleString('en')}
                </Accordion.Body>
            </Accordion.Item>
            <Accordion.Item eventKey={'1'}>
                <Accordion.Header><Badge>Indication</Badge></Accordion.Header>
                <Accordion.Body>
                    {props.drugInfo ? (
                        props.drugInfo.indication === 'N/A' ?
                            props.drugInfo.indication :
                            props.drugInfo.indication.split('.')[0] + '.'
                    ) : (
                        <Placeholder animation="glow">
                            <Placeholder xs={7} /> <Placeholder xs={4} /> <Placeholder xs={4} />{' '}
                            <Placeholder xs={6} /> <Placeholder xs={8} />
                        </Placeholder>
                    )}
                </Accordion.Body>
            </Accordion.Item>
            <Accordion.Item eventKey={'2'}>
                <Accordion.Header><Badge>Half-life</Badge></Accordion.Header>
                <Accordion.Body>
                    {props.drugInfo ? (
                        props.drugInfo.half_life
                    ) : (
                        <Placeholder animation="glow">
                            <Placeholder xs={7} /> <Placeholder xs={4} /> <Placeholder xs={4} />{' '}
                            <Placeholder xs={6} />
                        </Placeholder>
                    )}
                </Accordion.Body>
            </Accordion.Item>
            <Accordion.Item eventKey={'3'}>
                <Accordion.Header><Badge>Brands</Badge></Accordion.Header>
                <Accordion.Body>
                    {props.drugInfo ? (
                        // Format the brand list, so it is separated by commas and whitespaces
                        props.drugInfo.brands.map((brand, index) => {
                            return (
                                <span key={index}>
                                    {brand}
                                    {index !== props.drugInfo.brands.length - 1 ? ', ' : ''}
                                </span>
                            )
                        })
                    ) : (
                        <Placeholder animation="glow">
                            <Placeholder xs={7} /> <Placeholder xs={4} /> <Placeholder xs={4} />{' '}
                            <Placeholder xs={6} />
                        </Placeholder>
                    )}
                </Accordion.Body>
            </Accordion.Item>
        </Accordion>
    );
}

export default DrugAccordion;