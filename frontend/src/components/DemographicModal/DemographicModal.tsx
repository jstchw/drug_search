import { Modal, Col, Row, Placeholder } from "react-bootstrap";
import { useTermDataBatch } from "../../hooks/useTermDataBatch";
import { useUrlParams } from "../../hooks/useUrlParams";
import { URLParams } from "../../types";
import { searchAgeGroups, searchSex } from "../../constants";
import SimpleTermChart from "../ApexChart/SimpleTermChart";
import { Bug } from "@phosphor-icons/react";
import useDemographicStore from "../../stores/demographicStore";

const getParamsArray = (
  term: string,
  searchBy: string,
  searchMode: string,
): URLParams[] => {
  const paramList: URLParams[] = [];

  Object.entries(searchAgeGroups).forEach(([_, ageGroup]) => {
    searchSex.forEach((sexOption) => {
      const params: URLParams = {
        terms: [term],
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

const DemographicModal = () => {
  // Get the states from the store
  const [show, setShow] = useDemographicStore((state) => [
    state.showDemographic,
    state.setShowDemographic,
  ]);

  const term = useDemographicStore((state) => state.demographicTerm);

  const searchBy = useDemographicStore((state) => state.demographicType);

  const {
    params: { searchMode },
  } = useUrlParams();



  const paramsArray = getParamsArray(term, searchBy, searchMode);

  const { paramDataArray, error, isLoading } = useTermDataBatch(paramsArray);

  return (
    <Modal show={show} onHide={() => setShow(!show)} centered size={"xl"}>
      <Modal.Header closeButton>
        <Modal.Title>{term}</Modal.Title>
        <div className={"mx-2 vr"} />
        <span className={"text-secondary"}>Demographic breakdown</span>
      </Modal.Header>
      <Modal.Body>
        <Row>
          {paramDataArray !== undefined
            ? Object.entries(paramDataArray).map(([key, value]) => {
                return (
                  <Col key={key} md={6}>
                    <Row
                      className={
                        "d-flex align-items-center justify-content-center text-center"
                      }
                    >
                      <Row
                        className={"rounded shadow p-2"}
                        style={{ width: "fit-content" }}
                      >
                        <div className={"fs-4"}>{value.params["Sex"]}</div>
                        <div className={"fs-5"}>{value.params["Age"]}</div>
                      </Row>
                    </Row>
                    <Row>
                      <SimpleTermChart data={value.data} />
                    </Row>
                  </Col>
                );
              })
            : isLoading
              ? Array.from({ length: 6 }).map((_, i) => (
                  <Col key={i} md={6} className={"my-2"}>
                    <ChartPlaceholder />
                  </Col>
                ))
              : error && (
                  <Col md={12} className={"text-center"}>
                    <Bug
                      className={"mx-auto display-1 text-secondary"}
                      weight={"light"}
                    />
                    <div className={"display-4"}>Sorry!</div>
                    <div className={"display-6 text-secondary"}>
                      We couldn't find any data for this term.
                    </div>
                  </Col>
                )}
        </Row>
      </Modal.Body>
    </Modal>
  );
};

export default DemographicModal;
