import React from "react";
import {Col, Container, Row} from "react-bootstrap";
import PieChart from "../PieChart/PieChart";
import {CapsulePill, CapsulePill as CapsulePillIcon} from "react-bootstrap-icons";
import DrugDescription from "../DrugDescription";


const SearchResultObject = ({ searchResults, searchTerm }) => {
    const { termCountDict, totalCount } = searchResults.result

    console.log('SearchResultObject re-rendered', searchTerm)
    return (
        <div>
            <Container className="m-5">
                <Row>
                    <Col className="border border-dark rounded py-2">
                        <p className="text-lg-start" style={{ border: "1px solid black" }}><CapsulePillIcon />{searchTerm}</p>
                        <DrugDescription drugName={searchTerm} />
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