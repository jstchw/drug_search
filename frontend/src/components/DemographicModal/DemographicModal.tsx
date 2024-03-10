import { Modal, Col, Row, Nav, Button } from 'react-bootstrap';
import { searchAgeGroups, searchSex } from '../../constants';
import DemographicComparsionChart from '../ApexChart/DemographicComparsionChart';
import useDemographicStore from '../../stores/demographicStore';
import React from 'react';
import DonutChart from '../ApexChart/DonutChart';
import ToggleSwitch from '../ToggleSwitch/ToggleSwitch';
import { motion } from 'framer-motion';
import { mapParamToLabel } from '../../utils/utils';
import { searchTypes } from '../../constants';

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

  const demographicType = useDemographicStore((state) => state.demographicType);

  const [hasFdaDemographicData, setHasFdaDemographicData] = React.useState<boolean>(false);
  const handleFdaDemographicStatusChange = (status: boolean) => {
    setHasFdaDemographicData(status);
  };

  const [isDemographicChartAdvanced, setIsDemographicChartAdvanced] = React.useState<boolean>(false);
  const toggleDemographicChartAdvanced = React.useCallback(() => {
    setIsDemographicChartAdvanced((prev) => !prev);
  }, []);

  const [hasFdaDistributionData, setHasFdaDistributionData] = React.useState<boolean>(false);
  const handleFdaDistributionStatus = (status: boolean) => {
    setHasFdaDistributionData(status);
  };

  const [hasPmDistributionData, setHasPmDistributionData] = React.useState<boolean>(false);
  const handlePmDistributionStatus = (status: boolean) => {
    setHasPmDistributionData(status);
  };

  const [hasPmDemographicData, setHasPmDemographicData] = React.useState<boolean>(false);
  const handlePmDemographicStatusChange = (status: boolean) => {
    setHasPmDemographicData(status);
  };

  // Aggregatation states
  const [aggregateType, setAggregateType] = React.useState<AggregateType>('Sex'); // default

  const [currentPageKey, setCurrentPageKey] = React.useState<string>('');

  React.useEffect(() => {
    setCurrentPageKey(defaultPageKeys[aggregateType]);
  }, [aggregateType]);

  const [isPubmedDistribution, setIsPubmedDistribution] = React.useState<boolean>(false);
  const togglePubmedDistribution = React.useCallback(() => {
    setIsPubmedDistribution((prev) => !prev);
  }, []);

  const [demographicDataSource, setDemographicDataSource] = React.useState<'fda' | 'pm'>('fda');
  const toggleDemographicDataSource = React.useCallback(() => {
    setDemographicDataSource((prev) => (prev === 'fda' ? 'pm' : 'fda'));
  }, []);

  /* Framer motion variants */
  const containerVariants = {
    hidden: { opacity: 0 },
    visible: {
      opacity: 1,
      transition: {
        type: 'spring',
        stiffness: 260,
        damping: 20,
        staggerChildren: 0.1,
      },
    },
  };

  const itemVariants = {
    hidden: { opacity: 0 },
    visible: { opacity: 1 },
  };

  return (
    <Modal show={show} onHide={handleShow} centered size={'xl'}>
      <Modal.Header closeButton>
        <Button variant={'outline-primary me-2'} style={{ pointerEvents: 'none'}}>{mapParamToLabel(demographicType, searchTypes)}</Button>
        <Modal.Title>{term}</Modal.Title>
        <div className={'mx-2 vr'} />
        <span className={'text-secondary'}>Demographic breakdown</span>
      </Modal.Header>
      <Modal.Body>
        {hasFdaDistributionData && (
          <>
            <Row className={'mb-2 text-center'}>
              <span className={'fs-2 fw-light'}>Report statistics</span>
            </Row>
            <Row className={hasFdaDistributionData ? 'd-flex' : 'd-none'}>
              <ToggleSwitch handleToggleSwitch={togglePubmedDistribution}>View Pubmed data</ToggleSwitch>
            </Row>
          </>
        )}
        <motion.div
          key={'distribution'}
          layout
          variants={containerVariants}
          initial={'hidden'}
          animate={'visible'}
          exit={{ opacity: 0 }}
        >
          <Row>
            <Col className={'mb-3 text-center'}>
              <Row>
                <DonutChart source={'fda'} type={'age_group'} onDataStatusChange={handleFdaDistributionStatus} />
              </Row>
              {isPubmedDistribution && (
                <motion.div variants={itemVariants}>
                  <Row className={'mt-4'}>
                    <DonutChart source={'pm'} type={'age_group'} onDataStatusChange={handlePmDistributionStatus} />
                  </Row>
                </motion.div>
              )}
            </Col>
            <Col className={'mb-3 text-center'}>
              <Row>
                <DonutChart source={'fda'} type={'patient_sex'} onDataStatusChange={handleFdaDistributionStatus} />
              </Row>
              {isPubmedDistribution && (
                <motion.div variants={itemVariants}>
                  <Row className={'mt-4'}>
                    <DonutChart source={'pm'} type={'patient_sex'} onDataStatusChange={handlePmDistributionStatus} />
                  </Row>
                </motion.div>
              )}
            </Col>
          </Row>
          <motion.div variants={itemVariants}>
            {hasFdaDemographicData && (
              <div>
                <Row className={'text-center mt-3 mb-1'}>
                  <span className={'fs-2 fw-light'}>Demographic breakdown</span>
                </Row>
                <Row className={'mb-2'}>
                  <ToggleSwitch handleToggleSwitch={toggleDemographicChartAdvanced}>Advanced view</ToggleSwitch>
                </Row>
                {demographicType === 'side_effect' && (
                  <Row className={'mb-2'}>
                    <ToggleSwitch handleToggleSwitch={toggleDemographicDataSource}>View Pubmed data</ToggleSwitch>
                  </Row>
                )}
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
                onDataStatusChange={
                  demographicDataSource === 'fda' ? handleFdaDemographicStatusChange : handlePmDemographicStatusChange
                }
                advancedView={isDemographicChartAdvanced}
                source={demographicDataSource}
              />
            </Row>
          </motion.div>
        </motion.div>
      </Modal.Body>
    </Modal>
  );
};

export default DemographicModal;
