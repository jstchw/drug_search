import { Modal, Col, Row, Placeholder, Nav } from "react-bootstrap";
import { useTermDataBatch } from "../../hooks/useTermDataBatch";
import { useUrlParams } from "../../hooks/useUrlParams";
import { URLParams, DemographicData } from "../../types";
import { searchAgeGroups, searchSex } from "../../constants";
import SimpleTermChart from "../ApexChart/SimpleTermChart";
import { Bug } from "@phosphor-icons/react";
import useDemographicStore from "../../stores/demographicStore";
import { Carousel } from "react-responsive-carousel";
import React from "react";

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

const filterData = (data: DemographicData[], filterType: keyof DemographicData['params']) => {
  const aggreagateData: Record<string, DemographicData[]> = {};

  data.forEach((entry) => {
    const key = entry.params[filterType];
    if (key) {
      aggreagateData[key] = aggreagateData[key] || [];
      aggreagateData[key]!.push(entry);
    }
  });

  return aggreagateData;
}

const ChartPlaceholder = () => {
  return (
    <Placeholder className={"d-flex flex-column align-items-center"} as="div" animation="glow">
      <Placeholder as={"h2"} xs={4} bg="dark" />
      <Placeholder as={"div"} xs={12} bg="dark" style={{ height: "45vw" }} />
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

  const [aggregateType, setAggregateType] = React.useState<string>("Age");

  const aggregatedData = filterData(paramDataArray || [], aggregateType);

  const [carouselIndex, setCarouselIndex] = React.useState<number>(0);

  return (
    <Modal show={show} onHide={() => setShow(!show)} centered size={"xl"}>
      <Modal.Header closeButton>
        <Modal.Title>{term}</Modal.Title>
        <div className={"mx-2 vr"} />
        <span className={"text-secondary"}>Demographic breakdown</span>
        <div className={"mx-2 vr"} />
        <Nav
          variant={"pills"}
          className="justify-content-center align-items-center"
          onSelect={(index) => {
            setAggregateType(index as string)
            setCarouselIndex(0)
          }}
          defaultActiveKey={aggregateType}
        >
          <span className={"me-2"}>Group by:</span>
          <Nav.Item>
            <Nav.Link eventKey={"Age"}>Age</Nav.Link>
          </Nav.Item>
          <Nav.Item>
            <Nav.Link eventKey={"Sex"}>Sex</Nav.Link>
          </Nav.Item>
        </Nav>
      </Modal.Header>
      <Modal.Body>
        <Row>
          {Object.entries(aggregatedData).length !== 0 && (
            <>
            <Nav 
            variant={"pills"} 
            className="justify-content-center mb-3" 
            onSelect={(index) => setCarouselIndex(parseInt(index as string))}
            activeKey={carouselIndex}
            >
              {Object.keys(aggregatedData).map((key, index) => (
                <Nav.Item key={index}>
                  <Nav.Link eventKey={index}>{key}</Nav.Link>
                </Nav.Item>
              ))}
            </Nav>
              <Carousel
                showThumbs={false}
                showIndicators={false}
                showArrows={false}
                showStatus={false}
                selectedItem={carouselIndex}
                swipeable={false}
              >
                {Object.values(aggregatedData).map((data, index) => (
                  <Row key={index}>
                    <SimpleTermChart data={data} aggregateType={aggregateType}/>
                  </Row>
                ))}
              </Carousel>
            </>
            )
          }
          { isLoading
              ? 
              <ChartPlaceholder />
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
