import { Modal, Col, Row, Nav } from 'react-bootstrap';
import { searchAgeGroups, searchSex } from '../../constants';
import DemographicComparsionChart from '../ApexChart/DemographicComparsionChart';
import useDemographicStore from '../../stores/demographicStore';
import React from 'react';
import DonutChart from '../ApexChart/DonutChart';

type AggregateType = 'Sex' | 'Age';

const defaultPageKeys: Record<AggregateType, string> = {
  Sex: searchSex[0]!.label,
  Age: Object.keys(searchAgeGroups)[0]!,
};

const DemographicModal = () => {
  // Get the states from the store
  const [show, setShow] = useDemographicStore((state) => [state.showDemographic, state.setShowDemographic]);

  const term = useDemographicStore((state) => state.demographicTerm);

  const groupPageKeys = useDemographicStore((state) => state.groupPageKeys);

  const [hasDemographicData, setHasDemographicData] = React.useState<boolean>(false);

  const handleDemographicStatusChange = (status: boolean) => {
    setHasDemographicData(status);
  };

  const [hasDistributionData, setHasDistributionData] = React.useState<boolean>(false);

  const handleDistributionStatusChange = (status: boolean) => {
    setHasDistributionData(status);
  };

  // Aggregatation states
  const [aggregateType, setAggregateType] = React.useState<AggregateType>('Sex'); // default

  const [currentPageKey, setCurrentPageKey] = React.useState<string>('');

  React.useEffect(() => {
    setCurrentPageKey(defaultPageKeys[aggregateType]);
  }, [aggregateType]);

  return (
    <Modal show={show} onHide={() => setShow(!show)} centered size={'xl'}>
      <Modal.Header closeButton>
        <Modal.Title>{term}</Modal.Title>
        <div className={'mx-2 vr'} />
        <span className={'text-secondary'}>Demographic breakdown</span>
      </Modal.Header>
      <Modal.Body>
        {hasDistributionData && (
          <Row className={'mb-3 text-center'}>
            <span className={'fs-2 fw-light'}>Report statistics</span>
          </Row>
        )}
        <Row>
          <Col className={'mb-3 text-center'}>
            <DonutChart type={'age_group'} onDataStatusChange={handleDistributionStatusChange} />
          </Col>
          <Col className={'mb-3 text-center'}>
            <DonutChart type={'patient_sex'} onDataStatusChange={handleDistributionStatusChange} />
          </Col>
        </Row>
        {hasDemographicData && (
          <div>
            <Row className={'text-center my-4'}>
              <span className={'fs-2 fw-light'}>Demographic breakdown</span>
            </Row>
            <Row className={'text-center mb-2'}>
              <span className={'text-secondary'}>Group by:</span>
            </Row>
            <Row>
              <Nav
                variant={'pills'}
                className="justify-content-center align-items-center mb-3"
                onSelect={(index) => {
                  setAggregateType(index as AggregateType);
                }}
                defaultActiveKey={aggregateType}
              >
                <Nav.Item className={'mx-1'}>
                  <Nav.Link eventKey={'Sex'}>Sex</Nav.Link>
                </Nav.Item>
                <Nav.Item className={'mx-1'}>
                  <Nav.Link eventKey={'Age'}>Age</Nav.Link>
                </Nav.Item>
              </Nav>
            </Row>
            <Row>
              <Nav variant={'pills'} className="justify-content-center mb-3">
                {groupPageKeys.map((key, index) => (
                  <Nav.Item key={index} className={'mx-1'}>
                    <Nav.Link
                      eventKey={index}
                      active={currentPageKey === key}
                      onClick={() => setCurrentPageKey(key)}
                      className={'outline-primary'}
                    >
                      {key}
                    </Nav.Link>
                  </Nav.Item>
                ))}
              </Nav>
            </Row>
          </div>
        )}
        <Row>
          <DemographicComparsionChart
            aggregateType={aggregateType}
            currentPageKey={currentPageKey}
            onDataStatusChange={handleDemographicStatusChange}
          />
        </Row>
      </Modal.Body>
    </Modal>
  );
};

export default DemographicModal;
