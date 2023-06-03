import React from "react";
import {Col, Container, Row} from "react-bootstrap";
import PieChart from "../PieChart/PieChart";
import './SearchResultObject.css'
import DrugDescription from "../DrugDescription/DrugDescription";


const SearchResultObject = ({ searchResults, searchTerm }) => {
    const { termCountDict, totalCount } = searchResults.result
    const [retrievedTerm, setRetrievedTerm] = React.useState(null)

    return (
        <div>
            <Container className="m-4">
                {searchResults && <h1>{retrievedTerm}</h1>}
            </Container>
            <Container className="m-5">
                <Row>
                    <Col className="drug-desc-col">
                        <DrugDescription drugName={searchTerm} onRetrieved={setRetrievedTerm}/>
                    </Col>
                    <Col>
                        <p>Adverse effects for all groups</p>
                        {searchResults && <PieChart termCountDict={termCountDict} totalCount={totalCount} />}
                    </Col>
                </Row>
            </Container>
        </div>
    );
}

export default SearchResultObject;