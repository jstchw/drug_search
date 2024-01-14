import { Card, Badge } from 'react-bootstrap';
import './DrugDescription.css';
import {DrugProperties} from "../../types";


const DrugDescription = (props: {drugInfo: DrugProperties}) => {
    return (
        <Card className={'drug-description-card mb-4'}>
            {props.drugInfo.molecule_url &&
                <Card.Img src={props.drugInfo.molecule_url}
                          alt={props.drugInfo.name}
                          className='molecule'/>}

            <Card.Body className="text-start">
                {props.drugInfo && props.drugInfo.iupac &&
                    <p><Badge>IUPAC</Badge>&nbsp;<span className={'iupac'}>{props.drugInfo.iupac}</span></p>
                }
                {props.drugInfo && props.drugInfo.classification &&
                    <p><Badge>CLASS</Badge>&nbsp;{props.drugInfo.classification}</p>
                }
                {props.drugInfo && props.drugInfo.formula &&
                    <p><Badge>FORMULA</Badge>&nbsp;{props.drugInfo.formula}</p>
                }
            </Card.Body>
        </Card>
    );
};

export default DrugDescription;
