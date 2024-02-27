import { Modal, Col, Row, Nav } from 'react-bootstrap';
import { searchAgeGroups, searchSex } from '../../constants';
import DemographicComparsionChart from '../ApexChart/DemographicComparsionChart';
import useDemographicStore from '../../stores/demographicStore';
import React from 'react';
import DonutChart from '../ApexChart/DonutChart';
import PubmedSwitch from '../ChartSection/PubmedSwitch';
import { motion } from 'framer-motion';

type AggregateType = 'Sex' | 'Age';

const defaultPageKeys: Record<AggregateType, string> = {
  Sex: searchSex[0]!.label,
  Age: Object.keys(searchAgeGroups)[0]!,
};

const DemographicModal = () => {
  // Get the states from the store
  const [show, setShow] = useDemographicStore((state) => [state.showDemographic, state.setShowDemographic]);

  const handleShow = () => {
    setShow(!show);
  };

  const term = useDemographicStore((state) => state.demographicTerm);

  const groupPageKeys = useDemographicStore((state) => state.groupPageKeys);

  const [hasDemographicData, setHasDemographicData] = React.useState<boolean>(false);

  const handleDemographicStatusChange = (status: boolean) => {
    setHasDemographicData(status);
  };

  const [hasFdaDistributionData, setHasFdaDistributionData] = React.useState<boolean>(false);
  const handleFdaDistributionStatus = (status: boolean) => {
    setHasFdaDistributionData(status);
  };

  const [hasPmDistributionData, setHasPmDistributionData] = React.useState<boolean>(false);
  const handlePmDistributionStatus = (status: boolean) => {
    setHasPmDistributionData(status);
  };

  // Aggregatation states
  const [aggregateType, setAggregateType] = React.useState<AggregateType>('Sex'); // default

  const [currentPageKey, setCurrentPageKey] = React.useState<string>('');

  React.useEffect(() => {
    setCurrentPageKey(defaultPageKeys[aggregateType]);
  }, [aggregateType]);

  const [isPubmedAgeDistribution, setIsPubmedAgeDistribution] = React.useState<boolean>(false);
  const togglePubmedAgeDistribution = React.useCallback(() => {
    setIsPubmedAgeDistribution((prev) => !prev);
  }, []);

  const [isPubmedSexDistribution, setIsPubmedSexDistribution] = React.useState<boolean>(false);
  const togglePubmedSexDistribution = React.useCallback(() => {
    setIsPubmedSexDistribution((prev) => !prev);
  }, []);

  return (
    <Modal show={show} onHide={handleShow} centered size={'xl'}>
      <Modal.Header closeButton>
        <Modal.Title>{term}</Modal.Title>
        <div className={'mx-2 vr'} />
        <span className={'text-secondary'}>Demographic breakdown</span>
      </Modal.Header>
      <Modal.Body>
          {hasFdaDistributionData && (
            <Row className={'mb-3 text-center'}>
              <span className={'fs-2 fw-light'}>Report statistics</span>
            </Row>
          )}
          <motion.div
            key={'distribution'}
            layout
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            transition={{
              type: 'spring',
              stiffness: 260,
              damping: 20,
            }}
            exit={{ opacity: 0 }}
            >
          <Row>
            <Col className={'mb-3 text-center'}>
              <Row className={hasFdaDistributionData ? 'd-flex' : 'd-none'}>
                <PubmedSwitch handlePubmedSwitch={togglePubmedAgeDistribution} />
              </Row>
              <DonutChart source={'fda'} type={'age_group'} onDataStatusChange={handleFdaDistributionStatus} />
              {isPubmedAgeDistribution && (
                <DonutChart source={'pm'} type={'age_group'} onDataStatusChange={handlePmDistributionStatus} />
              )}
            </Col>
            <Col className={'mb-3 text-center'}>
              <Row className={hasFdaDistributionData ? 'd-flex' : 'd-none'}>
                <PubmedSwitch handlePubmedSwitch={togglePubmedSexDistribution} />
              </Row>
              <DonutChart source={'fda'} type={'patient_sex'} onDataStatusChange={handleFdaDistributionStatus} />
              {isPubmedSexDistribution && (
                <DonutChart source={'pm'} type={'patient_sex'} onDataStatusChange={handlePmDistributionStatus} />
              )}
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
                      <Nav.Item key={`${key}-${index}`} className={'mx-1'}>
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
          </motion.div>
      </Modal.Body>
    </Modal>
  );
};

export default DemographicModal;
