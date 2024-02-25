import { Accordion } from 'react-bootstrap';
import { ReadMore } from '../ReadMore/ReadMore';
import { Key } from 'react';
import { DrugProperties } from '../../types';
import { PencilSimple, ClockCountdown, Pill } from '@phosphor-icons/react';

const DrugAccordion = (props: { drugInfo: DrugProperties }) => {
  return (
    <Accordion>
      {props.drugInfo && props.drugInfo.indication && (
        <Accordion.Item eventKey={'0'}>
          <Accordion.Header>
            <div className={'d-flex justify-content-center fs-5'}>
              <PencilSimple weight={'light'} />
              <div className={'vr mx-2'} />
              <span>Indication</span>
            </div>
          </Accordion.Header>
          <Accordion.Body>
            {props.drugInfo.indication === 'N/A'
              ? props.drugInfo.indication
              : props.drugInfo.indication.split('.')[0] + '.'}
          </Accordion.Body>
        </Accordion.Item>
      )}
      {props.drugInfo && props.drugInfo.half_life && (
        <Accordion.Item eventKey={'1'}>
          <Accordion.Header>
            <div className={'d-flex justify-content-center fs-5'}>
              <ClockCountdown weight={'light'} />
              <div className={'vr mx-2'} />
              <span>Half-Life</span>
            </div>
          </Accordion.Header>
          <Accordion.Body>{props.drugInfo.half_life}</Accordion.Body>
        </Accordion.Item>
      )}
      {props.drugInfo && props.drugInfo.brands && (
        <Accordion.Item eventKey={'2'}>
          <Accordion.Header>
            <div className={'d-flex justify-content-center fs-5'}>
              <Pill weight={'light'} />
              <div className={'vr mx-2'} />
              <span>Brands</span>
            </div>
          </Accordion.Header>
          <Accordion.Body>
            <ReadMore>
              {props.drugInfo.brands.map((brand: string, index: Key) => (
                <span key={index}>
                  {brand}
                  {index !== (props.drugInfo.brands?.length ?? 0) - 1 ? ', ' : ' '}
                </span>
              ))}
            </ReadMore>
          </Accordion.Body>
        </Accordion.Item>
      )}
    </Accordion>
  );
};

export default DrugAccordion;
