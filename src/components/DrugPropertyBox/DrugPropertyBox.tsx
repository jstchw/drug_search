import {DrugProperties} from "../../types";
import DrugGroups from "../DrugGroups/DrugGroups";
import {Col, Row} from "react-bootstrap";
import DrugDescription from "../DrugDescription/DrugDescription";
import DrugAccordion from "../DrugAccordion/DrugAccordion";

const DrugPropertyBox = (props: {drug: DrugProperties, isSingle: boolean}) => {

    console.log(props.drug)
    return (
        <>
            <Col style={{ textAlign: 'center' }}>
                <h1>{props.drug.name}</h1>
                <Row className="justify-content-center">
                    <Col xs={12} md={8}>
                        <DrugGroups drugGroups={props.drug.groups || []}/>
                    </Col>
                </Row>
                {props.isSingle ? (
                    // For a single drug, put DrugDescription and DrugAccordion side by side
                    <Row className="justify-content-center">
                        <Col>
                            <DrugDescription drugInfo={props.drug}/>
                        </Col>
                        <Col>
                            <DrugAccordion drugInfo={props.drug}/>
                        </Col>
                    </Row>
                ) : (
                    // For multiple drugs, stack DrugDescription and DrugAccordion
                    <>
                        <Row className="justify-content-center">
                            <Col xs={12} md={8}>
                                <DrugDescription drugInfo={props.drug}/>
                            </Col>
                        </Row>
                        <Row className="justify-content-center">
                            <Col xs={12} md={8}>
                                <div style={{maxWidth: '100%', overflow: 'auto'}}>
                                    <DrugAccordion drugInfo={props.drug}/>
                                </div>
                            </Col>
                        </Row>
                    </>
                )}
            </Col>
        </>


    )
}

export default DrugPropertyBox;