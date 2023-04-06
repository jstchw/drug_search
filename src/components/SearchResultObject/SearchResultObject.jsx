import React from "react";
import {Col, Container, Row} from "react-bootstrap";
import PieChart from "../PieChart/PieChart";


const SearchResultObject = ({ searchResults, searchTerm }) => {
    const { termCountDict, totalCount } = searchResults.result

    return (
        <div>
            <h2>{searchTerm}</h2>
            <Container>
                <Row>
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