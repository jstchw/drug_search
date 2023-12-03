import {Accordion, Badge} from "react-bootstrap";
import {ReadMore} from "../ReadMore/ReadMore";
import {Key} from "react";
import { DrugInfo } from "../../types";

const DrugAccordion = (props: {drugInfo: DrugInfo}) => {
    return (
        <Accordion>
            {props.drugInfo && props.drugInfo.indication && (
                <Accordion.Item eventKey={'0'}>
                    <Accordion.Header><Badge>Indication</Badge></Accordion.Header>
                    <Accordion.Body>
                        {props.drugInfo.indication === 'N/A' ?
                            props.drugInfo.indication :
                            props.drugInfo.indication.split('.')[0] + '.'}
                    </Accordion.Body>
                </Accordion.Item>
            )}
            {props.drugInfo && props.drugInfo.half_life && (
                <Accordion.Item eventKey={'1'}>
                    <Accordion.Header><Badge>Half-life</Badge></Accordion.Header>
                    <Accordion.Body>
                        {props.drugInfo.half_life}
                    </Accordion.Body>
                </Accordion.Item>
            )}
            {props.drugInfo && props.drugInfo.brands && (
                <Accordion.Item eventKey={'2'}>
                    <Accordion.Header><Badge>Brands</Badge></Accordion.Header>
                    <Accordion.Body>
                        <ReadMore>
                            {props.drugInfo.brands.map((brand: string, index: Key) => (
                                <span key={index}>{brand}{index !== props.drugInfo.brands.length - 1 ? ', ' : ' '}</span>
                            ))}
                        </ReadMore>
                    </Accordion.Body>
                </Accordion.Item>
            )}
        </Accordion>
    );
}

export default DrugAccordion;