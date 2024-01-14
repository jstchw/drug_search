import {DrugProperties} from "../../types";
import DrugGroups from "../DrugGroups/DrugGroups";
import {Col, Row} from "react-bootstrap";
import DrugDescription from "../DrugDescription/DrugDescription";
import DrugAccordion from "../DrugAccordion/DrugAccordion";

const DrugPropertyBox = (props: {drug: DrugProperties}) => {

    return (
        <>
            <Col style={{ textAlign: 'center' }}>
                <h1>{props.drug.name}</h1>
                <Row className="justify-content-center">
                    <Col xs={12} md={8}>
                        <DrugGroups drugGroups={props.drug.groups || []}/>
                    </Col>
                </Row>
                <Row className="justify-content-center">
                    <Col xs={12} md={8}>
                        <DrugDescription drugInfo={props.drug} />
                    </Col>
                </Row>
                <Row className="justify-content-center">
                    <Col xs={12} md={8}>
                        <div style={{ maxWidth: '100%', overflow: 'auto' }}>
                            <DrugAccordion drugInfo={props.drug} />
                        </div>
                    </Col>
                </Row>
            </Col>
        </>

    )
}

export default DrugPropertyBox;