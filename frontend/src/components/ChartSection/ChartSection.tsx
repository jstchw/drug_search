import React from "react";
import CasesTimeSeriesChart from "../ApexChart/CasesTimeSeriesChart";
import TermCarousel from "../ApexChart/TermCarousel";
import "react-responsive-carousel/lib/styles/carousel.min.css";
import { Row, Col, ToggleButton } from "react-bootstrap";
import { Clock } from "@phosphor-icons/react";
import { isMobile } from "react-device-detect";
import { useAreParamsFiltered } from "../../hooks/useAreParamsFiltered";
import { getChartWarning } from "../ApexChart/TermCarousel";
import { useUrlParams } from "../../hooks/useUrlParams";

type ViewUnfilteredButtonProps = {
  noFilterRequest: boolean;
  setNoFilterRequest: React.Dispatch<React.SetStateAction<boolean>>;
  id: string;
  value: string;
};

const viewUnfilteredButton = ({
  noFilterRequest,
  setNoFilterRequest,
  id,
  value,
}: ViewUnfilteredButtonProps) => (
  <div className={"d-flex justify-content-center align-items-center"}>
    <ToggleButton
      id={id}
      value={value}
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
  const { params } = useUrlParams();

  const [noFilterTimeSeriesRequest, setNoFilterTimeSeriesRequest] =
    React.useState<boolean>(false);

  const [noFilterTermCarouselRequest, setNoFilterTermCarouselRequest] =
    React.useState<boolean>(false);

  const [childrenRendered, setChildrenRendered] = React.useState({
    timeSeries: false,
    termCarousel: false,
  });

  const handleChildRender = (child: string) => {
    setChildrenRendered((prev) => ({ ...prev, [child]: true }));
  };

  return (
    <div className={"mt-4"}>
      {/* Time series section */}
      <Row className={"d-flex justify-content-center"}>
        <Col
          xs={isMobile || noFilterTimeSeriesRequest ? 12 : 6}
          className="mb-4"
        >
          {/* If the children are rendered, show the title */}
          {childrenRendered.timeSeries && (
            <>
              <Row>
                <h3
                  className={"d-flex justify-content-center align-items-center"}
                >
                  <Clock weight={"light"} className={"text-secondary"} />
                  <div className={"vr mx-2"} />
                  Reports over time
                </h3>
              </Row>

              {areParamsFiltered &&
                viewUnfilteredButton({
                  noFilterRequest: noFilterTimeSeriesRequest,
                  setNoFilterRequest: setNoFilterTimeSeriesRequest,
                  id: "unfilteredTimeSeries",
                  value: "unfilteredTimeSeries",
                })}
            </>
          )}

          <Row className={"mb-4"}>
            <Col>
              {noFilterTimeSeriesRequest && (
                <h5 className={"d-flex justify-content-center mt-3"}>
                  Filtered data
                </h5>
              )}
              <CasesTimeSeriesChart
                onRender={() => handleChildRender("timeSeries")}
              />
            </Col>

            {noFilterTimeSeriesRequest && (
              <Col>
                <h5 className={"d-flex justify-content-center mt-3"}>
                  Unfiltered data
                </h5>
                <CasesTimeSeriesChart
                  noFilterRequest={true}
                  onRender={() => handleChildRender("timeSeries")}
                />
              </Col>
            )}
          </Row>
        </Col>
      </Row>

      {/* Term carousel section */}
      <Row className={"d-flex justify-content-center mb-4"}>
        <Col
          xs={isMobile || noFilterTermCarouselRequest ? 12 : 6}
          className="mb-4"
        >
          {/* If the children are rendered, show the title */}
          {childrenRendered.termCarousel && (
            <>
              <Row>
                <h3 className={"d-flex justify-content-center"}>
                  {getChartWarning(params)}
                </h3>
              </Row>

              {areParamsFiltered &&
                viewUnfilteredButton({
                  noFilterRequest: noFilterTermCarouselRequest,
                  setNoFilterRequest: setNoFilterTermCarouselRequest,
                  id: "unfilteredTermCarousel",
                  value: "unfilteredTermCarousel",
                })}
            </>
          )}

          <Row>
            <Col xs={noFilterTermCarouselRequest ? 6 : 12}>
              {noFilterTermCarouselRequest && (
                <h5 className={"d-flex justify-content-center mt-3"}>
                  Filtered data
                </h5>
              )}
              <TermCarousel
                onRender={() => handleChildRender("termCarousel")}
              />
            </Col>

            {noFilterTermCarouselRequest && (
              <Col xs={6}>
                <h5 className={"d-flex justify-content-center mt-3"}>
                  Unfiltered data
                </h5>
                <TermCarousel
                  noFilterRequest={true}
                  onRender={() => handleChildRender("termCarousel")}
                />
              </Col>
            )}
          </Row>
        </Col>
      </Row>
    </div>
  );
};

export default ChartSection;
