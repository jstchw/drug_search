import { Modal, Col, Row, Placeholder } from "react-bootstrap";
import { useTermDataBatch } from "../../hooks/useTermDataBatch";
import { useUrlParams } from "../../hooks/useUrlParams";
import { URLParams } from "../../types";
import { searchAgeGroups, searchSex } from "../../constants";
import SimpleTermChart from "../ApexChart/SimpleTermChart";
import React from "react";

const getParamsArray = (
  term: string[],
  searchBy: string,
  searchMode: string,
): URLParams[] => {
  const paramList: URLParams[] = [];

  Object.entries(searchAgeGroups).forEach(([_, ageGroup]) => {
    searchSex.forEach((sexOption) => {
      const params: URLParams = {
        terms: term,
        searchBy,
        searchMode,
        sex: sexOption.param,
        age: {
          min: ageGroup.min.toString(),
          max: ageGroup.max.toString(),
        },
      };
      paramList.push(params);
    });
  });

  return paramList;
};

const ChartPlaceholder = () => {
  return (
    <Placeholder as="div" animation="glow">
      <Placeholder as={"h2"} xs={4} bg="dark" />
      <Placeholder as={"div"} xs={12} bg="dark" style={{ height: "20vw" }} />
    </Placeholder>
  );
};

interface DemographicModalProps {
  show: boolean;
  setShow: (show: boolean) => void;
  term: string[];
}

const DemographicModal: React.FC<DemographicModalProps> = ({
  show,
  setShow,
  term,
}) => {
  const {
    params: { searchBy, searchMode },
  } = useUrlParams();

  const paramsArray = getParamsArray(term, searchBy, searchMode);

  const { paramDataArray, error, loading } = useTermDataBatch(paramsArray);

  return (
    <Modal show={show} onHide={() => setShow(!show)} centered size={"xl"}>
      <Modal.Header closeButton>
        <Modal.Title>{term[0]}</Modal.Title>
        <div className={"mx-2 vr"} />
        <span className={"text-secondary"}>Demographic breakdown</span>
      </Modal.Header>
      <Modal.Body>
        <Row>
          {paramDataArray !== null
            ? Object.entries(paramDataArray).map(([key, value]) => {
                return (
                  <Col key={key} md={6}>
                    <SimpleTermChart paramDataValue={value} />
                  </Col>
                );
              })
            : loading &&
              Array.from({ length: 6 }).map((_, i) => (
                <Col key={i} md={6} className={"my-2"}>
                  <ChartPlaceholder />
                </Col>
              ))}
        </Row>
      </Modal.Body>
    </Modal>
  );
};

export default DemographicModal;
