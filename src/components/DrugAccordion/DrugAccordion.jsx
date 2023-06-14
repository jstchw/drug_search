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
                    {props.retrievedTermArr ? (
                        props.retrievedTermArr.indication === 'N/A' ?
                            props.retrievedTermArr.indication :
                            props.retrievedTermArr.indication.split('.')[0] + '.'
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
                    {props.retrievedTermArr ? (
                        props.retrievedTermArr.half_life
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
                    {props.retrievedTermArr ? (
                        // Format the brand list, so it is separated by commas and whitespaces
                        props.retrievedTermArr.brands.map((brand, index) => {
                            return (
                                <span key={index}>
                                    {brand}
                                    {index !== props.retrievedTermArr.brands.length - 1 ? ', ' : ''}
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