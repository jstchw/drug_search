import React from 'react';
import CasesTimeSeriesChart from '../ApexChart/CasesTimeSeriesChart';
import TermCarousel from '../ApexChart/TermCarousel';
import 'react-responsive-carousel/lib/styles/carousel.min.css';
import { Row, Col, ToggleButton } from 'react-bootstrap';
import { Clock } from '@phosphor-icons/react';
import { isMobile } from 'react-device-detect';
import { useAreParamsFiltered } from '../../hooks/useAreParamsFiltered';
import { getChartWarning } from '../ApexChart/TermCarousel';
import { useUrlParams } from '../../hooks/useUrlParams';
import { AnimatePresence } from 'framer-motion';
import PubmedSwitch from './PubmedSwitch';

type ViewUnfilteredButtonProps = {
  noFilterRequest: boolean;
  setNoFilterRequest: React.Dispatch<React.SetStateAction<boolean>>;
  id: string;
  value: string;
};

const viewUnfilteredButton = ({ noFilterRequest, setNoFilterRequest, id, value }: ViewUnfilteredButtonProps) => (
  <div className={'d-flex justify-content-center align-items-center'}>
    <ToggleButton
      id={id}
      value={value}
      variant={'outline-primary'}
      className={'d-flex justify-content-center align-items-center my-2'}
      type={'checkbox'}
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

  const [noFilterTimeSeriesRequest, setNoFilterTimeSeriesRequest] = React.useState<boolean>(false);
  const [isPubmedTimeSeries, setIsPubmedTimeSeries] = React.useState<boolean>(false);
  const togglePubmedIncludedTimeSeries = React.useCallback(() => {
    setIsPubmedTimeSeries((prev) => !prev);
  }, []);

  const [noFilterTermCarouselRequest, setNoFilterTermCarouselRequest] = React.useState<boolean>(false);

  const [childrenRendered, setChildrenRendered] = React.useState({
    timeSeries: false,
    termCarousel: false,
  });

  const handleChildRender = (child: string) => {
    setChildrenRendered((prev) => ({ ...prev, [child]: true }));
  };

  return (
    <div className={'mt-4'}>
      {/* Time series section */}
      <Row className={'d-flex justify-content-center'}>
        <Col xs={isMobile || noFilterTimeSeriesRequest || isPubmedTimeSeries ? 12 : 6} className="mb-4">
          {/* If the children are rendered, show the title */}
          {childrenRendered.timeSeries && (
            <>
              <Row>
                <h3 className={'d-flex justify-content-center align-items-center'}>
                  <Clock weight={'light'} className={'text-secondary'} />
                  <div className={'vr mx-2'} />
                  Reports over time
                </h3>
              </Row>

              {/* If the params are not filtered, show the Pubmed switch */}
              {/* There is no point in showing the button when params are filtered becasue the data is too sparse */}
              {!areParamsFiltered && (
              <Row className={'my-1'}>
                <PubmedSwitch handlePubmedSwitch={togglePubmedIncludedTimeSeries}/>
              </Row>
              )}

              {areParamsFiltered &&
                viewUnfilteredButton({
                  noFilterRequest: noFilterTimeSeriesRequest,
                  setNoFilterRequest: setNoFilterTimeSeriesRequest,
                  id: 'unfilteredTimeSeries',
                  value: 'unfilteredTimeSeries',
                })}
            </>
          )}

          <AnimatePresence>
            <Row className={'mb-4'}>
              <Col>
                {noFilterTimeSeriesRequest && <h5 className={'d-flex justify-content-center mt-3'}>Filtered data</h5>}
                {isPubmedTimeSeries && <h5 className={'d-flex justify-content-center mt-3'}>FAERS data</h5>}
                <CasesTimeSeriesChart key={'genericTimeSeries'} onRender={() => handleChildRender('timeSeries')} />
              </Col>

              {isPubmedTimeSeries && (
                <Col>
                  <h5 className={'d-flex justify-content-center mt-3'}>Pubmed data</h5>
                  <CasesTimeSeriesChart
                    key={'pubmedTimeSeries'}
                    isPubmed={true}
                    onRender={() => handleChildRender('timeSeries')}
                  />
                </Col>
              )}

              {noFilterTimeSeriesRequest && (
                <Col>
                  <h5 className={'d-flex justify-content-center mt-3'}>Unfiltered data</h5>
                  <CasesTimeSeriesChart
                    key={'unfilteredTimeSeries'}
                    noFilterRequest={true}
                    onRender={() => handleChildRender('timeSeries')}
                  />
                </Col>
              )}
            </Row>
          </AnimatePresence>
        </Col>
      </Row>

      {/* Term carousel section */}
      <Row className={'d-flex justify-content-center mb-4'}>
        <Col xs={isMobile || noFilterTermCarouselRequest ? 12 : 6} className="mb-4">
          {/* If the children are rendered, show the title */}
          {childrenRendered.termCarousel && (
            <>
              <Row>
                <h3 className={'d-flex justify-content-center'}>{getChartWarning(params)}</h3>
              </Row>
              <Row className={'d-flex justify-content-center text-secondary'}>
                Percentage is calculated using only the presented data.
              </Row>

              {areParamsFiltered &&
                viewUnfilteredButton({
                  noFilterRequest: noFilterTermCarouselRequest,
                  setNoFilterRequest: setNoFilterTermCarouselRequest,
                  id: 'unfilteredTermCarousel',
                  value: 'unfilteredTermCarousel',
                })}
            </>
          )}

          <AnimatePresence>
            <Row>
              <Col xs={noFilterTermCarouselRequest ? 6 : 12}>
                {noFilterTermCarouselRequest && <h5 className={'d-flex justify-content-center mt-3'}>Filtered data</h5>}
                <TermCarousel onRender={() => handleChildRender('termCarousel')} />
              </Col>

              {noFilterTermCarouselRequest && (
                <Col xs={6}>
                  <h5 className={'d-flex justify-content-center mt-3'}>Unfiltered data</h5>
                  <TermCarousel noFilterRequest={true} onRender={() => handleChildRender('termCarousel')} />
                </Col>
              )}
            </Row>
          </AnimatePresence>
        </Col>
      </Row>
    </div>
  );
};

export default ChartSection;
