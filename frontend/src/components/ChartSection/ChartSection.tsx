import React from "react";
import CasesTimeSeriesChart from "../ApexChart/CasesTimeSeriesChart";
import TermCarousel from "../ApexChart/TermCarousel";
import "react-responsive-carousel/lib/styles/carousel.min.css";
import { Row, Col, ToggleButton } from "react-bootstrap";
import { Clock } from "@phosphor-icons/react";
import { isMobile } from "react-device-detect";
import { useAreParamsFiltered } from "../../hooks/useAreParamsFiltered";

type ViewUnfilteredButtonProps = {
  noFilterRequest: boolean;
  setNoFilterRequest: React.Dispatch<React.SetStateAction<boolean>>;
};

const viewUnfilteredButton = ({
  noFilterRequest,
  setNoFilterRequest,
}: ViewUnfilteredButtonProps) => (
  <div className={"d-flex justify-content-center align-items-center"}>
    <ToggleButton
      id={"noFilterTimeSeriesRequestButton"}
      value={"noFilterTimeSeriesRequest"}
      variant={"outline-primary"}
      className={"d-flex justify-content-center align-items-center"}
      type={"checkbox"}
      checked={noFilterRequest}
      onChange={() => setNoFilterRequest(!noFilterRequest)}
    >
      View unfiltered data
    </ToggleButton>
  </div>
);

const ChartSection = () => {
  const areParamsFiltered = useAreParamsFiltered();
  console.log(areParamsFiltered)

  const [noFilterTimeSeriesRequest, setNoFilterTimeSeriesRequest] =
    React.useState<boolean>(false);

  return (
    <div className={"mt-4"}>
      <Row className={"d-flex justify-content-center"}>
        <Col
          xs={isMobile || noFilterTimeSeriesRequest ? 12 : 6}
          className="mb-4"
        >
          <Row>
            <h3 className={"d-flex justify-content-center align-items-center"}>
              <Clock weight={"light"} className={"text-secondary"} />
              <div className={"vr mx-2"} />
              Reports over time
            </h3>
          </Row>
          {areParamsFiltered && viewUnfilteredButton({
            noFilterRequest: noFilterTimeSeriesRequest,
            setNoFilterRequest: setNoFilterTimeSeriesRequest,
          })}
          <Row className={"mb-4"}>
            <Col>
              {noFilterTimeSeriesRequest && (
                <h5 className={"d-flex justify-content-center mt-3"}>
                  Filtered data
                </h5>
              )}
              <CasesTimeSeriesChart />
            </Col>
            {noFilterTimeSeriesRequest && (
              <Col>
                <h5 className={"d-flex justify-content-center mt-3"}>
                  Unfiltered data
                </h5>
                <CasesTimeSeriesChart noFilterRequest={true} />
              </Col>
            )}
          </Row>
        </Col>
      </Row>
      <Row className={"d-flex justify-content-center mb-4"}>
        <Col xs={isMobile ? 12 : 6} className="mb-4">
          <TermCarousel />
        </Col>
      </Row>
    </div>
  );
};

export default ChartSection;
