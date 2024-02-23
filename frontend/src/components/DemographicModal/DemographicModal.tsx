import { Modal, Col, Row, Nav } from "react-bootstrap";
import { searchAgeGroups, searchSex } from "../../constants";
import SimpleTermChart from "../ApexChart/SimpleTermChart";
import useDemographicStore from "../../stores/demographicStore";
import React from "react";
import DonutChart from "../ApexChart/DonutChart";

type AggregateType = "Sex" | "Age";

const defaultPageKeys: Record<AggregateType, string> = {
  Sex: searchSex[0]!.label,
  Age: Object.keys(searchAgeGroups)[0]!,
};

const DemographicModal = () => {
  // Get the states from the store
  const [show, setShow] = useDemographicStore((state) => [
    state.showDemographic,
    state.setShowDemographic,
  ]);

  const term = useDemographicStore((state) => state.demographicTerm);

  const groupPageKeys = useDemographicStore((state) => state.groupPageKeys);

  // Aggregatation states
  const [aggregateType, setAggregateType] =
    React.useState<AggregateType>("Sex"); // default

  const [currentPageKey, setCurrentPageKey] = React.useState<string>("");

  React.useEffect(() => {
    setCurrentPageKey(defaultPageKeys[aggregateType]);
  }, [aggregateType]);

  return (
    <Modal show={show} onHide={() => setShow(!show)} centered size={"xl"}>
      <Modal.Header closeButton>
        <Modal.Title>{term}</Modal.Title>
        <div className={"mx-2 vr"} />
        <span className={"text-secondary"}>Demographic breakdown</span>
      </Modal.Header>
      <Modal.Body>
        <Row className={"mb-3 text-center"}>
          <span className={"fs-2 fw-light"}>Report statistics</span>
        </Row>
        <Row>
          <Col className={"mb-3 text-center"}>
            <span className={"fs-4 fw-normal"}>Age distribution</span>
            <DonutChart type={"age_group"} />
          </Col>
          <Col className={"mb-3 text-center"}>
            <span className={"fs-4 fw-normal"}>Sex distribution</span>
            <DonutChart type={"patient_sex"} />
          </Col>
        </Row>
        <Row className={"text-center my-4"}>
          <span className={"fs-2 fw-light"}>Demographic breakdown</span>
        </Row>
          <Row className={"text-center mb-2"}>
            <span className={"text-secondary"}>Group by:</span>
          </Row>
          <Row>
            <Nav
              variant={"pills"}
              className="justify-content-center align-items-center mb-3"
              onSelect={(index) => {
                setAggregateType(index as AggregateType);
              }}
              defaultActiveKey={aggregateType}
            >
              <Nav.Item className={"mx-1"}>
                <Nav.Link eventKey={"Sex"}>Sex</Nav.Link>
              </Nav.Item>
              <Nav.Item className={"mx-1"}>
                <Nav.Link eventKey={"Age"}>Age</Nav.Link>
              </Nav.Item>
            </Nav>
          </Row>
        <Row>
          <Nav variant={"pills"} className="justify-content-center mb-3">
            {groupPageKeys.map((key, index) => (
              <Nav.Item key={index} className={"mx-1"}>
                <Nav.Link
                  eventKey={index}
                  active={currentPageKey === key}
                  onClick={() => setCurrentPageKey(key)}
                  className={"outline-primary"}
                >
                  {key}
                </Nav.Link>
              </Nav.Item>
            ))}
          </Nav>
        </Row>
        <Row>
          <SimpleTermChart
            aggregateType={aggregateType}
            currentPageKey={currentPageKey}
          />
        </Row>
      </Modal.Body>
    </Modal>
  );
};

export default DemographicModal;
