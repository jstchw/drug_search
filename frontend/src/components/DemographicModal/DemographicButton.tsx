import { Button } from 'react-bootstrap';
import { IdentificationCard } from '@phosphor-icons/react';
import useDemographicStore from '../../stores/demographicStore';
import React from 'react';
import { useUrlParams } from '../../hooks/useUrlParams';

const DemographicButton: React.FC<{ term: string }> = ({ term }) => {
  const {
    params: { searchBy },
  } = useUrlParams();

  const [showDemographic, setShowDemographic] = useDemographicStore((state) => [
    state.showDemographic,
    state.setShowDemographic,
  ]);

  const setDemographicTerm = useDemographicStore((state) => state.setDemographicTerm);
  const setDemographicType = useDemographicStore((state) => state.setDemographicType);

  const handleShowDemographic = () => {
    setShowDemographic(!showDemographic);
    setDemographicTerm(term);
    setDemographicType(searchBy);
  };

  return (
    <Button
      variant="outline-primary"
      className={'mb-4 d-flex align-items-center'}
      style={{ width: 'fit-content' }}
      onClick={handleShowDemographic}
    >
      <IdentificationCard weight={'light'} className={'fs-3'} />
      <div className={'vr mx-2'} />
      Demographic breakdown
    </Button>
  );
};

export default DemographicButton;
