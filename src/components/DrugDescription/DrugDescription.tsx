import { Card, Badge } from 'react-bootstrap';
import './DrugDescription.css';
import {DrugInfo} from "../../types";


const DrugDescription = (props: {drugInfo: DrugInfo}) => {
    return (
        <Card className={'drug-description-card mb-4'} style={{border: "none"}}>
            {props.drugInfo.molecule_url && <Card.Img src={props.drugInfo.molecule_url} alt={props.drugInfo.drug_name}
                                                className="molecule" style={{ backgroundColor: 'transparent'}} />}

            <Card.Body className="text-start">
                {props.drugInfo && props.drugInfo.iupac &&
                    <p><Badge>IUPAC</Badge>&nbsp;<span className={'iupac'}>{props.drugInfo.iupac}</span></p>
                }
                {props.drugInfo && props.drugInfo.drug_class &&
                    <p><Badge>Class</Badge>&nbsp;{props.drugInfo.drug_class}</p>
                }
                {props.drugInfo && props.drugInfo.formula &&
                    <p><Badge>Formula</Badge>&nbsp;{props.drugInfo.formula}</p>
                }
            </Card.Body>
        </Card>
    );
};

export default DrugDescription;
