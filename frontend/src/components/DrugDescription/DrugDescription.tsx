import { Card } from 'react-bootstrap';
import './DrugDescription.css';
import { DrugProperties } from '../../types';

const DrugDescription = (props: { drugInfo: DrugProperties }) => {
  return (
    <Card className={'drug-description-card mb-4'}>
      {props.drugInfo.molecule_url && (
        <Card.Img src={props.drugInfo.molecule_url} alt={props.drugInfo.name} className="molecule" />
      )}

      <Card.Body className="text-start">
        {props.drugInfo && props.drugInfo.iupac && (
          <div className={'d-flex align-items-center mb-3'}>
            <span className={'fs-5'}>IUPAC</span>
            <div className={'vr mx-2'} />
            <span className={'iupac'}>{props.drugInfo.iupac}</span>
          </div>
        )}
        {props.drugInfo && props.drugInfo.classification && (
          <div className={'d-flex align-items-center mb-3'}>
            <span className={'fs-5'}>Class</span>
            <div className={'vr mx-2'} />
            {props.drugInfo.classification}
          </div>
        )}
        {props.drugInfo && props.drugInfo.formula && (
          <div className={'d-flex align-items-center'}>
            <span className={'fs-5'}>Formula</span>
            <div className={'vr mx-2'} />
            {props.drugInfo.formula}
          </div>
        )}
      </Card.Body>
    </Card>
  );
};

export default DrugDescription;
